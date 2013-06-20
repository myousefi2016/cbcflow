#!/usr/bin/env python

from headflow import *
from headflow.dol import *
set_log_level(100)


# List schemes
spd = ParamDict(
    u_degree=1,
    solver_p=("cg", "amg"),
    #solver_p=("gmres", "amg"),
    )
schemes = [
    IPCS(spd),
    #SegregatedIPCS(spd),
    #SegregatedIPCS_Optimized(spd),
    ]

# List problems
from cylinder import FlowAroundACylinder
from drivencavity import DrivenCavity
from beltrami import Beltrami
from pipe import Pipe

ppd = ParamDict(
    #dt = 1e-3,
    #T  = 1e-3 * 100,
    num_periods=0.02,
    )
problems = [
    #FlowAroundACylinder(ppd),
    #DrivenCavity(ppd),
    #Beltrami,
    Pipe(ppd),
    ]

# Configure postprocessing
pppd = ParamDict(casedir="pipe_results")
postprocessor = NSPostProcessor(pppd)

ppfield_pd = ParamDict(
    saveparams=ParamDict(
        save=True,
        ),
    timeparams=ParamDict(
        step_frequency=10,
        )
    )
wss = WSS(params=ppfield_pd)
velocity = Velocity(params=ppfield_pd)
pressure = Pressure(params=ppfield_pd)

postprocessor.add_fields([wss, velocity, pressure])

# Configure solver
npd = ParamDict(
    plot_solution=False,
    check_mem_usage=True,
    )

# Loop over schemes and problems
for scheme in schemes:
    for problem in problems:
        solver = NSSolver(problem, scheme, postprocessor, params=npd)
        solver.solve()
