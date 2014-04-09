#!/usr/bin/env python

import sys, os, itertools
import unittest

from cbcflow import *

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

def extract_table(data, refinement_levels, dts, fieldname):
    em = {}
    return { (refinement_level,dt): data.get((refinement_level,dt),em).get(fieldname)
             for refinement_level in refinement_levels for dt in dts }

# TODO: Find some better table formatting utils I have lying around somewhere
def print_table(caption, table, refinement_levels, dts):
    maxchars = max(len(v) for v in table.values())
    widthfmt = "%%%ds" % maxchars
    dtfmt = "%.4e"
    sep = "  "

    colheader = "dt       " + sep.join(widthfmt % (dtfmt % dt) for dt in dts)
    lines = ["", caption, colheader]

    for refinement_level in refinement_levels:
        row = "ref %5d  " % refinement_level
        row += sep.join(widthfmt % table[(refinement_level, dt)] for dt in dts)
        lines.append(row)

    table = "\n".join(lines)
    print table
    return table

class TestAnalyticalSolutionConvergence(DiscretizationSweepTestCase):

    def _make_fields(self):
        "Return postprocessing fields to apply in solve."
        # TODO: Make storing configurable, better in an automated test to have automatic in-memory analysis
        p = ParamDict(
            save=False,
            )
        fields = [
            Velocity(p),
            AnalyticalVelocity(p),
            ]

        xnorms = [
            DiffL2norm("Velocity", "AnalyticalVelocity", p),
            DiffH1norm("Velocity", "AnalyticalVelocity", p),
            DiffH1seminorm("Velocity", "AnalyticalVelocity", p),
            ]
        tnorms = []
        tnorms += [RunningMax(xn, p) for xn in xnorms]
        tnorms += [RunningL2norm(xn, p) for xn in xnorms]
        self._norm_field_names = [tn.name for tn in tnorms]

        fields += xnorms
        fields += tnorms
        return fields

    def _analyse_data(self, data):
        "Analyse the data provided by the discretization parameter sweep."
        # FIXME: Compare with reference data, and use assertions for values with known properties.
        #print data
        refinement_levels = self._refinement_levels()
        dts = self._dts()
        for fieldname in self._norm_field_names:
            table = extract_table(data, refinement_levels, dts, fieldname)
            table = map_dict_values(table, lambda v: "--" if v is None else ("%2.2e" % v))
            formatted = print_table(fieldname, table, refinement_levels, dts)

            # Store formatted table as a reference TODO: Store floats so we can use a precision argument for comparing
            self._write_reference(data["basedir"], fieldname, formatted)

            # While we're at it, check that we can read back the reference data
            readback = self._read_reference(data["basedir"], fieldname, basedir="output")
            self.assertEqual(readback.strip(), formatted.strip())


# FIXME: Need a better solution for importing demos
# Importing problems from cbcflow/demo/
# NB! Assuming run from the cbcflow/tests/ directory!
sys.path.insert(0, "../demo/undocumented/Poiseuille2D")
sys.path.insert(0, "../demo/undocumented/Poiseuille3D")
sys.path.insert(0, "../demo/undocumented/Womersley2D")
sys.path.insert(0, "../demo/undocumented/Womersley3D")
sys.path.insert(0, "../demo/undocumented/Beltrami")
from Poiseuille2D import Poiseuille2D
from Poiseuille3D import Poiseuille3D
from Womersley2D import Womersley2D
from Womersley3D import Womersley3D
from Beltrami import Beltrami


def load_tests(loader, standard_tests, none):

    tests = []

    # FIXME: Make fast and slow suite

    # FIXME: Add more schemes, use all of official_schemes
    schemes = [
    #    lambda: BottiPietro(),
    #    lambda: CoupledNonLinear(),
    #    lambda: CoupledPicard(),
    #    lambda: IPCS(),
    #    lambda: IPCS_Stabilized({'theta':0.0}),
        lambda: IPCS_Stabilized({'theta':0.5}),
    #    lambda: IPCS_Stabilized({'theta':1.0}),
        lambda: IPCS_Stable(),
    #    lambda: Karper(),
    #    lambda: PISO(),
    #    lambda: PenaltyIPCS(),
    #    lambda: SegregatedIPCS(),
    #    lambda: SegregatedIPCS_Optimized(),
    #    lambda: SegregatedPenaltyIPCS(),
    #    lambda: Stokes(),
        ]

    T = lambda dt: dt*5
    num_periods = None
    problems = [
        lambda refinement_level,dt: Poiseuille2D(ParamDict(refinement_level=refinement_level, dt=dt, T=None, num_periods=0.2)),
    #    lambda refinement_level,dt: Poiseuille3D(ParamDict(refinement_level=refinement_level, dt=dt, T=None, num_periods=0.2)),
        lambda refinement_level,dt: Womersley2D(ParamDict(refinement_level=refinement_level, dt=dt, T=None, num_periods=1.0)),
    #    lambda refinement_level,dt: Womersley3D(ParamDict(N=N, dt=dt, T=None, num_periods=1.0)), # FIXME: Parameterize by N
        lambda refinement_level,dt: Beltrami(ParamDict(refinement_level=refinement_level, dt=dt, T=1.0)),
        ]

    fast = True
    if fast:
        params = [dict(
            refinement_levels=[0,1],
            dts = [0.1, 0.05],
            )]
    elif 1:
        params = [dict(
            refinement_levels=[2,3],
            dts = [0.1, 0.05],
            )]
    else:
        params = [dict(
            refinement_levels=[3,4],
            dts = [0.01, 0.005, 0.0025, 0.00125],
            )]

    tests.append(make_suite(TestAnalyticalSolutionConvergence, [schemes, problems, params]))

    return unittest.TestSuite(tests)