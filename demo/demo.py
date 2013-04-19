#!/usr/bin/env python

# Hack to run without installing, useful while working
import sys
sys.path.insert(0,"../site-packages")

from dolfin import *
set_log_level(100)

from headflow import ParamDict, NSSolver

# Setup scheme
#from headflow import IPCS as Scheme
from headflow import SegregatedIPCS as Scheme
spd = ParamDict(
    dt=1e-4,
    u_degree=1,
    )
scheme = Scheme(spd)

# Setup problem
#from drivencavity import DrivenCavity as Problem
from beltrami import Beltrami as Problem
#from cylinder import FlowAroundACylinder as Problem
ppd = ParamDict(
    N=20,
    )
problem = Problem(ppd)

# Setup solver
npd = ParamDict(
    plot_solution=False,
    check_mem_usage=True,
    )
solver = NSSolver(problem, scheme, params=npd)

# Execute!
solver.solve()
