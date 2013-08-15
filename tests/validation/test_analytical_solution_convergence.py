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
    maxchars = max(len(v) for v in table.values())
    widthfmt = "%%%ds" % maxchars
    dtfmt = "%.4e"
    sep = "  "

    colheader = "dt       " + sep.join(widthfmt % (dtfmt % dt) for dt in dts)
    lines = ["", caption, colheader]

    for N in Ns:
        row = "N %5d  " % N
        row += sep.join(widthfmt % table[(N, dt)] for dt in dts)
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
        Ns = self._Ns()
        dts = self._dts()
        for fieldname in self._norm_field_names:
            table = extract_table(data, Ns, dts, fieldname)
            table = map_dict_values(table, lambda v: "--" if v is None else ("%2.2e" % v))
            formatted = print_table(fieldname, table, Ns, dts)
            self._write_reference(fieldname, formatted)

# Importing problems from headflow/demo/
sys.path.insert(0, "../demo") # NB! Assuming run from the headflow/tests/ directory!
#sys.path.insert(0, os.path.abspath( os.path.join(os.path.dirname(__file__), '../../demo/') ))
from pouseille2d import Pouseille2D
from pouseille3d import Pouseille3D
from womersley2d import Womersley2D
from womersley3d import Womersley3D
from beltrami import Beltrami

def load_tests(loader, standard_tests, none):

    tests = []

    # FIXME: Make fast and slow suite

    # FIXME: Add more schemes
    schemes = [
    #    lambda N,dt: BottiPietro(),
    #    lambda N,dt: CoupledNonLinear(),
    #    lambda N,dt: CoupledPicard(),
    #    lambda N,dt: IPCS(),
    #    lambda N,dt: IPCS_Stabilized(),
    #    lambda N,dt: IPCS_Stable({'theta':0.0}),
        lambda N,dt: IPCS_Stable({'theta':0.5}),
    #    lambda N,dt: IPCS_Stable({'theta':1.0}),
    #    lambda N,dt: Karper(),
    #    lambda N,dt: PISO(),
    #    lambda N,dt: PenaltyIPCS(),
    #    lambda N,dt: SegregatedIPCS(),
    #    lambda N,dt: SegregatedIPCS_Optimized(),
    #    lambda N,dt: SegregatedPenaltyIPCS(),
    #    lambda N,dt: Stokes(),
        ]

    problems = [ # FIXME: Limiting T for debugging:
        #lambda N,dt: Pouseille2D(ParamDict(N=N, dt=dt, T=dt*5, num_periods=None)),
        #lambda N,dt: Pouseille3D(ParamDict(N=N, dt=dt, T=dt*2, num_periods=None)),
        lambda N,dt: Womersley2D(ParamDict(N=N, dt=dt, T=None, num_periods=.1)),
        #lambda N,dt: Womersley3D(ParamDict(N=N, dt=dt, T=dt*2, num_periods=None)),
        #lambda N,dt: Beltrami(ParamDict(N=N, dt=dt, T=dt*2)),
        ]

    fast = False
    if fast:
        params = [dict(
            Ns  = [2, 4],
            dts = [0.1, 0.05],
            )]
    else:
        params = [dict(
        #    Ns  = [16, 32],
            Ns  = [8, 16],
            dts = [0.01, 0.005],
        #    dts = [0.01, 0.005, 0.0025, 0.00125],
            )]

    tests.append(make_suite(TestAnalyticalSolutionConvergence, [schemes, problems, params]))

    return unittest.TestSuite(tests)
