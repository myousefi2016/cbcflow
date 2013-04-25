#!/usr/bin/env python

# Hack to run without installing, useful while working
import sys
sys.path.insert(0,"../site-packages")

from headflow import NSSolver
import headflow.dol as dol

from channel import Problem
from headflow.schemes.ipcs_opt import Scheme

if __name__ == "__main__":

    # TODO: Parse commandline options to overwrite default params, use utils in ParamDict class

    # Set global DOLFIN parameters
    options = {
        "debug": False,
        "krylov_solver_absolute_tolerance": 1e-25,
        "krylov_solver_relative_tolerance": 1e-12,
        "krylov_solver_monitor_convergence": False,
        }
    dol.set_log_active(options["debug"])
    dol.parameters["form_compiler"]["cpp_optimize"] = True
    dol.parameters["krylov_solver"]["absolute_tolerance"] = options["krylov_solver_absolute_tolerance"]
    dol.parameters["krylov_solver"]["relative_tolerance"] = options["krylov_solver_relative_tolerance"]
    dol.parameters["krylov_solver"]["monitor_convergence"] = options["krylov_solver_monitor_convergence"]

    # FIXME: Restarting! Extract logic from ns script. Make part of NSSolver?

    # Create instance of problem
    p_params = Problem.default_params()
    problem = Problem(p_params)

    # Create instance of postprocessor
    #pp_params = PostProcessor.default_params()
    postproc = None #PostProcessor(pp_params)

    # Create instance of scheme
    s_params = Scheme.default_params()
    scheme = Scheme(s_params)

    # Instantiate generic solver
    ns_params = NSSolver.default_params()
    solver = NSSolver(problem,
                      scheme=scheme,
                      postprocessor=postproc,
                      params=ns_params)
    solver.solve()

