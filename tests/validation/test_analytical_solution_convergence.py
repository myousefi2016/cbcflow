#!/usr/bin/env python

import sys, os, itertools
import unittest

from headflow import *

from math import sqrt
import dolfin
dolfin.parameters["allow_extrapolation"] = True

from .discretization_sweep_test_case import DiscretizationSweepTestCase, make_suite


class TestAnalyticalSolutionConvergence(DiscretizationSweepTestCase):

    def _make_fields(self):
        "Return postprocessing fields to apply in solve."
        # FIXME: Make storing configurable, better in an automated test to have automatic in-memory analysis
        save = False
        p1 = ParamDict(
            save=True,
            )
        p3 = ParamDict(
            save=True,
            stride_timestep=10,
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

        print data

        # TODO: Move these utils to shared code?

        def map_dict_values(table, mapping):
            return { k:mapping(v) for k,v in table.iteritems() }

        def multi_get(top, *keys):
            em = {}
            d = top
            for k in keys[:-1]:
                d = d.get(k,em)
            return d.get(keys[-1])

        def extract_table(data, Ns, dts, fieldname, subfieldname):
            em = {}
            return { (N,dt): multi_get(data, (N,dt), fieldname, subfieldname)
                     for N in Ns for dt in dts }

        # TODO: Find some better table formatting utils I have lying around somewhere
        def print_table(table, Ns, dts):
            print
            print x
            print "dt    ",
            for dt in dts:
                print dt,
            for N in Ns:
                print "\nN  ", N,
                print '  '.join(table[(N, dt)] for dt in dts),
            print

        Ns = self._Ns()
        dts = self._dts()
        fieldname = "AnalyticalSolutionAnalyzer"
        subfields = ["u0", "u1", "u2", "p"]
        for x in subfields:
            table = extract_table(data, Ns, dts, fieldname, x)
            table = map_dict_values(table, lambda v: "--" if v is None else ("%2.2e" % v))
            print_table(table, Ns, dts)

# Importing problems from headflow/demo/
# NB! Assuming run from the headflow/tests/ directory!
sys.path.insert(0, "../demo")
from flow_around_cylinder import FlowAroundCylinder
from beltrami import Beltrami

def load_tests(loader, standard_tests, none):

    tests = []

    # FIXME: Make fast and slow suite

    # FIXME: Add more schemes
    schemes = [
        lambda: IPCS(),
        lambda: IPCS_Stable({'theta':0.5}),
        ]

    # FIXME: Add more problems
    problems = [
        lambda N,dt: Beltrami(ParamDict(N=N, dt=dt, T=dt*2)), # FIXME: Limiting T for debugging
        ]
    params = [dict(
    #    Ns  = [2, 4, 8, 16],
    #    dts = [0.1, 0.05, 0.025, 0.0125],
        Ns  = [2, 4],
        dts = [0.1, 0.05],
        )]
    tests.append(make_suite(TestAnalyticalSolutionConvergence, [schemes, problems, params]))

    return unittest.TestSuite(tests)
