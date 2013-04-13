# Hack to run without installing, useful while working
import sys
sys.path.insert(0,"../site-packages")

from dolfin import *
set_log_level(100)
from headflow import *

#from drivencavity import DrivenCavity
from beltrami import Beltrami as Problem
#from cylinder import FlowAroundACylinder as Problem
#from ipcs_segregated import Scheme


problem = Problem(ParamDict(N=20))

#Solver = IPCS
Solver = SegregatedIPCS

scheme = Solver(ParamDict(dt=1e-4, u_degree=1))

pd = ParamDict(
    plot_solution=False,
    check_mem_usage=True,
    )

solver = NSSolver(problem, scheme, params=pd)
solver.solve()





