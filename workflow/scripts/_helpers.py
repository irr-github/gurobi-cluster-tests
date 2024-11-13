import logging
import sys
from pathlib import Path
import os
import traceback

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


def configure_logging(snakemake, logger=None, skip_handlers=False, level="INFO"):
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
        logger.addHandler(logging.StreamHandler())
        logger.addHandler(logging.FileHandler(logfile))

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
