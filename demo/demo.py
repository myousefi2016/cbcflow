#!/usr/bin/env python

# Hack to run without installing, useful while working
import sys; sys.path.insert(0,"../site-packages")

from headflow import ParamDict, NSSolver
from headflow.dol import set_log_level
set_log_level(100)


# Configure scheme
#from headflow import IPCS
from headflow import SegregatedIPCS
spd = ParamDict(
    u_degree=1,
    )
scheme = SegregatedIPCS(spd)


# List problems
from cylinder import FlowAroundACylinder
from drivencavity import DrivenCavity
from beltrami import Beltrami
from pipe import Pipe
problems = [
    #FlowAroundACylinder,
    #DrivenCavity,
    #Beltrami,
    Pipe,
    ]

ppd = ParamDict(
    dt = 1e-4,
    T  = 1e-4 * 10,
    )

# Loop over problems
for Problem in problems:
    # Configure problem
    problem = Problem(ppd)

    # Setup solver
    npd = ParamDict(
        plot_solution=False,
        check_mem_usage=True,
        )
    solver = NSSolver(problem, scheme, params=npd)

    # Execute!
    solver.solve()

