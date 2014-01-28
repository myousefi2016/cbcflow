#!/usr/bin/env python

import sys, os, itertools
import unittest

from cbcflow import *

from math import sqrt
import dolfin
dolfin.parameters["allow_extrapolation"] = True

from discretization_sweep_test_case import DiscretizationSweepTestCase, make_suite


class TestQ(DiscretizationSweepTestCase):

    def _Ns(self):
        "Return range of spatial discretization parameters."
        return [8, 16]

    def _dts(self):
        "Return range of temporal discretization parameters."
        return [0.05, 0.025]

    def _make_fields(self):
        "Return postprocessing fields to apply in solve."
        return [Velocity(), Q(), Delta(), Lambda2()]

    def _analyse_data(self, data):
        "Analyse the data provided by the discretization parameter sweep."
        pass


# FIXME: Need a better solution for importing demos!
# Importing problems from cbcflow/demo/
# NB! Assuming run from the cbcflow/tests/ directory!
sys.path.insert(0, "../demo/undocumented/FlowAroundCylinder")
sys.path.insert(0, "../demo/undocumented/Beltrami")
from flow_around_cylinder import FlowAroundCylinder
from beltrami import Beltrami


def load_tests(loader, standard_tests, none):

    # FIXME: Make fast and slow suite

    # FIXME: Add more schemes, use all of official_schemes
    schemes = [
        lambda: IPCS(),
        lambda: IPCS_Stable({'theta':0.5}),
        ]

    # FIXME: Add more problems
    problems = [
        lambda N,dt: Beltrami(ParamDict(N=N, dt=dt, T=dt*2)), # FIXME: Limiting T for debugging
        lambda N,dt: FlowAroundCylinder(ParamDict(N=N, dt=dt, T=dt*2)), # FIXME: Limiting T for debugging
        ]

    return unittest.TestSuite([]) #make_suite(TestAnalyticalSolutionConvergence, [schemes, problems]) # FIXME: Enable when test case does something
