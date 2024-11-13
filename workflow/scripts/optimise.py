import os.path
import json
import yaml

# import sys,  traceback
import os

# order!
import logging
import gurobipy as grb

from os import PathLike, makedirs  # , environ
from _helpers import mock_snakemake, configure_logging

logger = logging.getLogger(__name__)


def load_gurobi_license(lic_path: PathLike) -> dict:
    """transform the WSL gurobi license file into a dictionary that can
    be used to set the gurobi env params. Not needed if env was set before

    Args:
        lic_path (PathLike): the path

    Returns:
        dict: a dictionary
    """
    # license looks like: ARG1=VAL1 ARG2=VAL2
    lic = yaml.load(open("/p/projects/rd3mod/gurobi_rc/gurobi.lic", "r"), Loader=yaml.FullLoader)
    return dict([tuple(x.split("=")) for x in lic.split(" ")])


def optimize(cfg_path: PathLike, solution_path: PathLike, gurobi_options: dict = None):
    """a toy optimisation model to test gurobi + snakemake on the cluster

    Args:
        cfg_path (PathLike): constraints for the problem (json)
        solution_path (PathLike): output path (json)
        gurobi_options (dict, optional): gurobi env options
    """
    with open(cfg_path) as f:
        data = json.load(f)

    if gurobi_options is None:
        gurobi_options = {}
    elif "LICENSEID" in gurobi_options:
        gurobi_options["LICENSEID"] = int(gurobi_options["LICENSEID"])

    logger.info(f"optimize Using Gurobi configuration: {cfg_path}")

    with grb.Env(params=gurobi_options) as env, grb.Model("simple_lp", env=env) as m:
        logger.info("Created Gurobi model")
        x = m.addVar(name="x", vtype=grb.GRB.CONTINUOUS)
        y = m.addVar(name="y", vtype=grb.GRB.CONTINUOUS)

        m.setObjective(x + y, grb.GRB.MAXIMIZE)

        # constraints
        m.addConstr(x + 2 * y <= data["constraint_1"], "c0")

        m.optimize()

        if m.status == grb.GRB.OPTIMAL:
            solution = {"x": x.x, "y": y.x, "objective": m.ObjVal}

    if not os.path.exists(os.path.dirname(solution_path)):
        makedirs(os.path.dirname(solution_path))

    with open(solution_path, "w") as f:
        json.dump(solution, f, indent=4)

    logger.info(f"Wrote solution to {solution_path}")


def main(license_path: PathLike, cfg_path: PathLike, solution_path: PathLike, threads=1):
    """A function that can be called by a snakemake rule 'run' instead of 'script'

    Args:
        license_path (PathLike): the path to the license (to set a gurobi env)
        cfg_path (PathLike): the opt problem config (constraint, json)
        solution_path (PathLike): output path
        threads (int, optional): n threads to solve the problem. Defaults to 1.
    """
    gurobi_cfg = load_gurobi_license(license_path)
    gurobi_cfg["Threads"] = threads

    logger.info(f"MAIN Using Gurobi configuration: {gurobi_cfg}")
    optimize(
        os.path.abspath(cfg_path),
        os.path.abspath(solution_path),
        gurobi_cfg,
    )


if __name__ == "__main__":
    # ===== this will be called by snakemake rule script: optimise.py =====

    # for debugging can emulate the snakemake script object created by snakemake
    # rule needs to match
    if "snakemake" not in globals():
        snakemake = mock_snakemake("solve")

    logger.info("PRE CONFIG TEST")
    configure_logging(snakemake, logger=logger, skip_handlers=False)
    logger.info("Reading license file")
    gurobi_cfg = load_gurobi_license(snakemake.input.license_file)
    gurobi_cfg["Threads"] = 1
    logger.info("Starting optimization")
    optimize(
        snakemake.input.constraints,
        snakemake.output[0],
        gurobi_cfg,
    )
    # optimize(snakemake.input.constraints, snakemake.output, gurobi_cfg)
