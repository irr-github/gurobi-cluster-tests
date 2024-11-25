import logging
import sys
from pathlib import Path
import os, sys
import traceback, subprocess

logger = logging.getLogger(__name__)

DEFAULT_TUNNEL_PORT = 1080


def setup_gurobi_tunnel_and_env(tunnel_config: dict, logger: logging.Logger = None):
    """A utility function to set up the Gurobi environment variables and establish an SSH tunnel on HPCs
    Otherwise the license check will fail if the compute nodes do not have internet access or a token server isn't set up

    Args:
        config (dict): the snakemake pypsa-china configuration
        logger (logging.Logger, optional): Logger. Defaults to None.
    """
    if tunnel_config.get("use_tunnel", False) is False:
        return
    logger.info("setting up tunnel")
    user = os.getenv("USER")  # User is pulled from the environment
    port = tunnel_config.get("port", DEFAULT_TUNNEL_PORT)
    ssh_command = f"ssh -fN -D {port} {user}@login01"

    try:
        # Run SSH in the background to establish the tunnel
        subprocess.Popen(ssh_command, shell=True)
        logger.info(f"SSH tunnel established on port {port}")
    # TODO don't handle unless neeeded
    except Exception as e:
        logger.error(f"Error starting SSH tunnel: {e}")
        sys.exit(1)

    os.environ["https_proxy"] = f"socks5://127.0.0.1:{port}"
    os.environ["SSL_CERT_FILE"] = "/p/projects/rd3mod/ssl/ca-bundle.pem_2022-02-08"
    os.environ["GRB_CAFILE"] = "/p/projects/rd3mod/ssl/ca-bundle.pem_2022-02-08"

    # Set up Gurobi environment variables
    os.environ["GUROBI_HOME"] = "/p/projects/rd3mod/gurobi1103/linux64"
    os.environ["PATH"] += f":{os.environ['GUROBI_HOME']}/bin"
    os.environ["LD_LIBRARY_PATH"] += f":{os.environ['GUROBI_HOME']}/lib"
    os.environ["GRB_LICENSE_FILE"] = "/p/projects/rd3mod/gurobi_rc/gurobi.lic"
    os.environ["GRB_CURLVERBOSE"] = "1"

    logger.info("Gurobi Environment variables & tunnel set up successfully.")


# from pypsa-eur code


def mock_snakemake(
    rulename: str,
    configfiles: list | str = None,
    **wildcards,
):
    """A function to enable scripts to run as standalone, giving them access to
     the snakefile rule input, outputs etc

    Args:
        rulename (str): the name of the rule
        configfiles (list or str, optional): the config file or config file list. Defaults to None.
        wildcards (optional):  keyword arguments fixing the wildcards (if any needed)

    Raises:

        FileNotFoundError: Config file not found

    Returns:
        snakemake.script.Snakemake: an object storing all the rule inputs/outputs etc
    """

    import snakemake as sm
    from snakemake.api import Workflow
    from snakemake.common import SNAKEFILE_CHOICES
    from snakemake.script import Snakemake
    from snakemake.settings.types import (
        ConfigSettings,
        DAGSettings,
        ResourceSettings,
        StorageSettings,
        WorkflowSettings,
    )

    for p in SNAKEFILE_CHOICES:
        if os.path.exists(p):
            snakefile = p
            break
    if configfiles is None:
        configfiles = []
    elif isinstance(configfiles, str):
        configfiles = [configfiles]

    resource_settings = ResourceSettings()
    config_settings = ConfigSettings(configfiles=map(Path, configfiles))
    workflow_settings = WorkflowSettings()
    storage_settings = StorageSettings()
    dag_settings = DAGSettings(rerun_triggers=[])
    workflow = Workflow(
        config_settings,
        resource_settings,
        workflow_settings,
        storage_settings,
        dag_settings,
        storage_provider_settings=dict(),
    )
    workflow.include(snakefile)

    if configfiles:
        for f in configfiles:
            if not os.path.exists(f):
                raise FileNotFoundError(f"Config file {f} does not exist.")
            workflow.configfile(f)

    workflow.global_resources = {}
    rule = workflow.get_rule(rulename)
    dag = sm.dag.DAG(workflow, rules=[rule])
    wc = wildcards
    job = sm.jobs.Job(rule, dag, wc)

    def make_accessable(*ios):
        for io in ios:
            for i, _ in enumerate(io):
                io[i] = os.path.abspath(io[i])

    make_accessable(job.input, job.output, job.log)
    snakemake = Snakemake(
        job.input,
        job.output,
        job.params,
        job.wildcards,
        job.threads,
        job.resources,
        job.log,
        job.dag.workflow.config,
        job.rule.name,
        None,
    )
    # create log and output dir if not existent
    for path in list(snakemake.log) + list(snakemake.output):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    return snakemake


def configure_logging(snakemake, logger=None, skip_handlers=False, level="INFO", fname=None):
    """
    Configure the basic behaviour for the logging module.
    Note: Must only be called once from the __main__ section of a script.
    The setup includes printing log messages to STDERR and to a log file defined
    by either (in priority order): snakemake.log.python, snakemake.log[0] or "logs/{rulename}.log".
    Additional keywords from logging.basicConfig are accepted via the snakemake configuration
    file under snakemake.config.logging.
    Parameters
    ----------
    snakemake : snakemake object
        Your snakemake object containing a snakemake.config and snakemake.log.
    skip_handlers : True | False (default)
        Do (not) skip the default handlers created for redirecting output to STDERR and file.
    """

    if not logger:
        logger = logging.getLogger()
        logger.info("Configuring logging")

    kwargs = snakemake.config.get("logging", dict())
    kwargs.setdefault("level", level)

    if skip_handlers is False:
        fallback_path = Path(__file__).parent.joinpath("..", "logs", f"{snakemake.rule}.log")
        default_logfile = snakemake.log[0] if snakemake.log else fallback_path
        logfile = snakemake.log.get("python", default_logfile)
        logger.setLevel(kwargs["level"])

        formatter = logging.Formatter("%(asctime)s - %(filename)s - %(levelname)s - %(message)s")
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

        file_handler = logging.FileHandler(logfile)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # def handle_exception(exc_type, exc_value, exc_traceback):
    #     # Log the exception
    #     logger = logging.getLogger()
    #     logger.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
    #     sys.excepthook = handle_exception

    def handle_exception(exc_type, exc_value, exc_traceback):
        # do not overload if KeyboardInterrupt
        if issubclass(exc_type, KeyboardInterrupt):
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
            return

        logger.error(
            "".join(
                [
                    "Uncaught exception: ",
                    *traceback.format_exception(exc_type, exc_value, exc_traceback),
                ]
            )
        )

    # Install exception handler
    sys.excepthook = handle_exception
