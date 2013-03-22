
import sys
sys.path.insert(0,"../site-packages")


import dolfin
from channel import Problem

if __name__ == "__main__":
    import sys
    from headflow import parse_cmdline_params
    options = parse_cmdline_params(sys.argv[1:])
    options["casedir"] = "tempcase"
    for key, value in options.iteritems():
        if key.startswith("solver.") and isinstance(value, str):
            options[key] = value.split(',')

    # Set global DOLFIN parameters
    dolfin.parameters["form_compiler"]["cpp_optimize"] = True
    dolfin.parameters["krylov_solver"]["absolute_tolerance"] = options["krylov_solver_absolute_tolerance"]
    dolfin.parameters["krylov_solver"]["relative_tolerance"] = options["krylov_solver_relative_tolerance"]
    dolfin.parameters["krylov_solver"]["monitor_convergence"] = options["krylov_solver_monitor_convergence"]

    # Set debug level
    dolfin.set_log_active(options["debug"])

    # Create instance of problem defined above
    problem = Problem(options)

    # For now just fetch the solver and start it manually:
    from headflow.solvers.ipcs_opt import Solver
    scheme = Solver(options)
    scheme.solve(problem)

    # ...

    # Later we should have a generic interface:
    from headflow import NSSolver
    solver = NSSolver(problem, params=options)

    solver.solve()
