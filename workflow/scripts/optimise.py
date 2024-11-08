import gurobipy as grb
import json

import os.path
from os import PathLike, makedirs

def optimize(cfg_path:PathLike, solution_path:PathLike):

    with open(cfg_path) as f:
        data = json.load(f)

    m = grb.Model("simple_lp")

    x  = m.addVar(name="x", vtype = grb.GRB.CONTINUOUS)
    y = m.addVar(name = "y", vtype=grb.GRB.CONTINUOUS)

    m.set_objective(x + y, grb.GRB.MAXIMIZE)

    # constraints
    m.set_constraint(x + 2*y <= data["constraint_1"], "c0")


    m.optimize()

    solution = {"x": x.X, "y": y.Y, "objective": m.objVal}

    if not os.path.exists(os.path.dirname(solution_path)):
        makedirs(os.path.dirname(solution_path))
    with open(solution_path, "w") as f:
        json.dump(  solution, f, indent=4)


if __name__ == "__main__":
    optimize(os.path.abspath("./data/constraints.json"), os.path.abspath("./data/solution.json"))
    # optimize(snakemake.input.constraints, snakemake.output)