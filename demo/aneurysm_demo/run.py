#!/usr/bin/env python

import sys; sys.path.insert(0,"../../site-packages")

from headflow import *
from headflow.dol import *
set_log_level(100)


# --- Configure problem
from aneurysm import DogAneurysm
problem_pd = ParamDict(
    num_periods=0.1,
    )
problem = DogAneurysm(problem_pd)


# --- Configure scheme
scheme_pd = ParamDict(
    u_degree=1,
    )
scheme = IPCS(scheme_pd)
#scheme = SegregatedIPCS(scheme_pd)
#scheme = SegregatedIPCS_Optimized(scheme_pd)


# -- Configure postprocessor
class PostProcessor(NSPostProcessor):
    def __init__(self, params=None):
        NSPostProcessor.__init__(self, params)

postprocessor = PostProcessor()

ppfield_pd = ParamDict(
    saveparams=ParamDict(
        save=True,
        ),
    timeparams=ParamDict(
        step_frequency=10,
        )
    )
#wss = WSS(params=ppfield_pd)
velocity = Velocity(params=ppfield_pd)
#pressure = Pressure(params=ppfield_pd)

postprocessor.add_field(velocity)
#postprocessor.add_fields([wss, velocity, pressure])

# --- Configure solver and run
pd = ParamDict(
    plot_solution=True,
    )
solver = NSSolver(problem, scheme, postprocessor, params=pd)
solver.solve()
