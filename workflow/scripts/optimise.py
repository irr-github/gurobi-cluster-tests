import gurobipy as grb
import json

import os.path
from os import PathLike, makedirs
import yaml

from snakemake import logging as snakemake_logging
import logging

# Get the Snakemake logger
logger = logging.getLogger("snakemake")


def load_gurobi_license(lic_path: PathLike) -> dict:
    """transform the WSL gurobi license file into a dictionary

    Args:
        lic_path (PathLike): the path

    Returns:
        dict: a dictionary
    """
    # license looks like: ARG1=VAL1 ARG2=VAL2
    lic = yaml.load(open("/p/projects/rd3mod/gurobi_rc/gurobi.lic", "r"), Loader=yaml.FullLoader)
    return dict([tuple(x.split("=")) for x in lic.split(" ")])


def optimize(cfg_path: PathLike, solution_path: PathLike, gurobi_options: dict):
    with open(cfg_path) as f:
        data = json.load(f)

    print("starting")

    gurobi_options["LICENSEID"] = int(gurobi_options["LICENSEID"])
    logger.info(f"optimize Using Gurobi configuration: {gurobi_cfg}")

    with grb.Env(params=gurobi_options) as env, grb.Model("simple_lp", env=env) as m:
        x = m.addVar(name="x", vtype=grb.GRB.CONTINUOUS)
        y = m.addVar(name="y", vtype=grb.GRB.CONTINUOUS)

        m.set_objective(x + y, grb.GRB.MAXIMIZE)

        # constraints
        m.set_constraint(x + 2 * y <= data["constraint_1"], "c0")

        m.optimize()

        solution = {"x": x.X, "y": y.Y, "objective": m.objVal}

    if not os.path.exists(os.path.dirname(solution_path)):
        makedirs(os.path.dirname(solution_path))
    with open(solution_path, "w") as f:
        json.dump(solution, f, indent=4)


def main(license_path: PathLike, cfg_path: PathLike, solution_path: PathLike, threads=1):
    gurobi_cfg = load_gurobi_license(license_path)
    gurobi_cfg["Threads"] = threads
    logger.info(f"MAIN Using Gurobi configuration: {gurobi_cfg}")
    optimize(
        os.path.abspath(cfg_path),
        os.path.abspath(solution_path),
        gurobi_cfg,
    )


if __name__ == "__main__":
    gurobi_cfg = {"Threads": 1}
    gurobi_cfg.update(load_gurobi_license(snakemake.input.license_file))
    # optimize(
    #     os.path.abspath("./data/constraints.json"),
    #     os.path.abspath("./data/solution.json"),
    #     gurobi_cfg,
    # )
    optimize(snakemake.input.constraints, snakemake.output, gurobi_cfg)
