#!/usr/bin/env python

import sys, os, itertools
import unittest

from headflow import *

from math import sqrt
import dolfin
dolfin.parameters["allow_extrapolation"] = True

from discretization_sweep_test_case import DiscretizationSweepTestCase, make_suite


# TODO: Move these utils to shared code?

def map_dict_values(container, mapping):
    "Construct a new dict with mapping function applied to each value in input dict."
    return { k: mapping(v) for k,v in container.iteritems() }

def multi_level_get(top, *keys):
    "Like top.get(key), but handling multi-level nested acces with given keys."
    em = {}
    d = top
    for k in keys[:-1]:
        d = d.get(k, em)
    return d.get(keys[-1])

def extract_table(data, Ns, dts, fieldname):
    em = {}
    return { (N,dt): data.get((N,dt),em).get(fieldname)
             for N in Ns for dt in dts }

# TODO: Find some better table formatting utils I have lying around somewhere
def print_table(caption, table, Ns, dts):
    print
    print caption
    print "dt    ",
    for dt in dts:
        print dt,
    for N in Ns:
        print "\nN  ", N,
        print '  '.join(table[(N, dt)] for dt in dts),
    print


class TestAnalyticalSolutionConvergence(DiscretizationSweepTestCase):

    def _make_fields(self):
        "Return postprocessing fields to apply in solve."
        # FIXME: Make storing configurable, better in an automated test to have automatic in-memory analysis
        save = False
        p1 = ParamDict(
            save=save,
            )
        p3 = ParamDict(
            save=save,
            stride_timestep=10,
            )
        return [
            #VelocityError(params=p1),
            #PressureError(params=p1),
            DiffL2norm("Velocity", "AnalyticalVelocity", params=p1),
            DiffH1norm("Velocity", "AnalyticalVelocity", params=p1),
            DiffH1seminorm("Velocity", "AnalyticalVelocity", params=p1),
            ]

    def _analyse_data(self, data):
        "Analyse the data provided by the discretization parameter sweep."
        # FIXME: Don't use print! Log to file and compare with reference data, and use assertions for values with known properties.
        self._print_table(data)

    def _print_table(self, data):
        print data
        Ns = self._Ns()
        dts = self._dts()
        fieldnames = [
            "DiffL2norm_Velocity_AnalyticalVelocity",
            "DiffH1norm_Velocity_AnalyticalVelocity",
            "DiffH1seminorm_Velocity_AnalyticalVelocity",
            ]
        #("L2norm_VelocityError", "H1norm_VelocityError", "H1seminorm_VelocityError", "L2norm_PressureError"):
        for fieldname in fieldnames:
            table = extract_table(data, Ns, dts, fieldname)
            table = map_dict_values(table, lambda v: "--" if v is None else ("%2.2e" % v))
            print_table(fieldname, table, Ns, dts)

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

    # FIXME: Add more problems: Pouseille2d, Pouseille3d, Womersley2d, Womersley3d
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
