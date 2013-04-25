#!/usr/bin/env python

# Hack to run without installing, useful while working
import sys; sys.path.insert(0,"../site-packages")

from headflow import *
from headflow.dol import *
set_log_level(100)


# List schemes
spd = ParamDict(
    u_degree=1,
    )
schemes = [
    SegregatedIPCS_Optimized(spd),
    ]

# List problems
from cylinder import FlowAroundACylinder
from drivencavity import DrivenCavity
from beltrami import Beltrami
from pipe import Pipe

ppd = ParamDict(
    #dt = 1e-3,
    #T  = 1e-3 * 100,
    )
problems = [
    #FlowAroundACylinder(ppd),
    #DrivenCavity(ppd),
    #Beltrami,
    Pipe(ppd),
    ]

# Configure solver
npd = ParamDict(
    plot_solution=False,
    check_mem_usage=True,
    )

# Loop over schemes and problems
for scheme in schemes:
    for problem in problems:
        solver = NSSolver(problem, scheme, params=npd)
        solver.solve()

