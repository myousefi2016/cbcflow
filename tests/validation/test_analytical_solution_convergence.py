#!/usr/bin/env python

import sys, os, itertools
import unittest

from headflow import *

from math import sqrt
import dolfin
dolfin.parameters["allow_extrapolation"] = True

from .discretization_sweep_test_case import DiscretizationSweepTestCase, make_suite


class TestAnalyticalSolutionConvergence(DiscretizationSweepTestCase):

    def _Ns(self):
        "Return range of spatial discretization parameters."
        return [2, 4, 8, 16]

    def _dts(self):
        "Return range of temporal discretization parameters."
        return [0.1, 0.05, 0.025, 0.0125]

    def _make_fields(self):
        "Return postprocessing fields to apply in solve."
        # FIXME: Make storing configurable, better in an automated test to have automatic in-memory analysis
        save = False
        p1 = ParamDict(
            saveparams=ParamDict(
                save=save,
                ),
            )
        p3 = ParamDict(
            saveparams=ParamDict(
                save=save,
                ),
            timeparams=ParamDict(
                step_frequency=10,
                )
            )
        return [
            AnalyticalSolutionAnalyzer(params=p1),
            EnergyAnalyzer(params=p1),
            WSS(params=p3),
            ]

    def _analyse_data(self, data):
        "Analyse the data provided by the discretization parameter sweep."
        # FIXME: Don't use print! Log to file and compare with reference data, and use assertions for values with known properties.
        self._print_table(data)

    def _print_table(self, data):
        # TODO: Split into computing and presenting table

        Ns = self._Ns()
        dts = self._dts()

        fieldname = "AnalyticalSolutionAnalyzer"
        subfields = ["u0", "u1", "u2", "p"]

        for x in subfields:
            print
            print
            print x
            print "dt ",
            for dt in dts:
                print dt,
            for N in Ns:
                print "\nN ", N,
                for dt in dts:
                    print " %2.2e " % data[(N, dt)][fieldname]["data"][x],


# Importing problems from headflow/demo/
# NB! Assuming run from the headflow/tests/ directory!
sys.path.insert(0, "../demo")
from flow_around_cylinder import FlowAroundCylinder
from beltrami import Beltrami

def load_tests(loader, standard_tests, none):

    # FIXME: Make fast and slow suite

    # FIXME: Add more schemes
    schemes = [
        lambda: IPCS(),
        lambda: IPCS_Stable({'theta':0.5}),
        ]

    # FIXME: Add more problems
    problems = [
        lambda N,dt: Beltrami(ParamDict(N=N, dt=dt, T=dt*2)), # FIXME: Limiting T for debugging
        lambda N,dt: FlowAroundCylinder(ParamDict(N=N, dt=dt, T=dt*2)), # FIXME: Limiting T for debugging
        ]

    return make_suite(TestAnalyticalSolutionConvergence, [schemes, problems])
