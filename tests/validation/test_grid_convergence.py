#!/usr/bin/env python

import sys, os, itertools
import unittest

from headflow import *

from math import sqrt
import dolfin
dolfin.parameters["allow_extrapolation"] = True

from .discretization_sweep_test_case import DiscretizationSweepTestCase, make_suite


class TestGridConvergence(DiscretizationSweepTestCase):

    def _Ns(self):
        "Return range of spatial discretization parameters."
        return [8, 16]

    def _dts(self):
        "Return range of temporal discretization parameters."
        return [0.05, 0.025]

    def _make_fields(self):
        "Return postprocessing fields to apply in solve."
        return [Velocity()]

    def _analyse_data(self, data):
        "Analyse the data provided by the discretization parameter sweep."
        # FIXME: Don't use print! Log to file and compare with reference data, and use assertions for values with known properties.
        self._analyze_spatial_convergence(data)
        self._analyze_temporal_convergence(data)

    def _analyze_spatial_convergence(self, data):
        "With time discretization fixed, loop over spatial discretization and compute norms."
        eps = 1e-14

        Ns = self._Ns()
        dts = self._dts()

        analysis = {}
        for dt in dts:
            analysis[dt] = []
            for N0, N1 in zip(Ns[:-1], Ns[1:]):
                data0 = data[(N0, dt)]
                data1 = data[(N1, dt)]
                a = { "dt": dt, "N0": N0, "N1": N1 }

                try:
                    u_coarse = data0["Velocity"]["data"]
                    u_fine = data1["Velocity"]["data"]
                    uc = dolfin.interpolate(u_coarse, u_fine.function_space())

                    a["u_coarse_norm"] = dolfin.assemble(u_coarse**2*dolfin.dx())
                    a["u_fine_norm"] = dolfin.assemble(u_fine**2*dolfin.dx())
                    a["u_diff_norm"] = dolfin.assemble((u_fine - uc)**2*dolfin.dx())

                except Exception as e:
                    a["exception"] = str(e)

                analysis[dt].append(a)

        # Slight modification of the prints that Kent did before:
        for dt in dts:
            for a in analysis[dt]:
                print "Difference between level ", a["N1"], " and ", a["N0"], " with dt ", dt, \
                  " is ", a["u_diff_norm"], " versus u_fine_norm ", a["u_fine_norm"], \
                  " relative root ", sqrt(a["u_diff_norm"] / (a["u_fine_norm"]+eps))

        # FIXME: Add assertions to validate convergence rate from analysis[dt][:]
        for dt in dts:
            rates = []
            for a0, a1 in zip(analysis[dt][:-1], analysis[dt][1:]):
                try:
                    # FIXME: Better rate estimation formula!
                    # Hack to avoid division by zero
                    e0 = sqrt(a0["u_diff_norm"] / (a0["u_fine_norm"]+eps))
                    e1 = sqrt(a1["u_diff_norm"] / (a1["u_fine_norm"]+eps))
                    rate = e0 / (e1+eps)
                except:
                    rate = None
                rates.append(rate)
            print "Convergence rates w.r.t N at dt=%g: %s" % (dt, rates)

    def _analyze_temporal_convergence(self, data):
        "With spatial discretization fixed, loop over time discretization and compute norms."
        eps = 1e-14

        Ns = self._Ns()
        dts = self._dts()

        analysis = {}
        for N in Ns:
            analysis[N] = []
            for dt0, dt1 in zip(dts[:-1], dts[1:]):
                data0 = data[(N, dt0)]
                data1 = data[(N, dt1)]
                a = { "dt0": dt0, "dt1": dt1, "N": N }

                try:
                    u_coarse = data0["Velocity"]["data"]
                    u_fine = data1["Velocity"]["data"]
                    uc = dolfin.interpolate(u_coarse, u_fine.function_space())

                    a["u_coarse_norm"] = dolfin.assemble(u_coarse**2*dolfin.dx())
                    a["u_fine_norm"] = dolfin.assemble(u_fine**2*dolfin.dx())
                    a["u_diff_norm"] = dolfin.assemble((u_fine - uc)**2*dolfin.dx())

                except Exception as e:
                    a["exception"] = str(e)

        # Slight modification of the prints that Kent did before:
        for N in Ns:
            for a in analysis[N]:
                print "Difference between dt ", a["dt1"], " and ", a["dt0"], " at level ", N, \
                  " is ", a["u_diff_norm"], " versus u_fine_norm ", a["u_fine_norm"], \
                  " relative root ", sqrt(a["u_diff_norm"] / (a["u_fine_norm"]+eps))

        # FIXME: Add assertions to validate convergence rate from analysis[N][:]
        for N in Ns:
            rates = []
            for a0, a1 in zip(analysis[N][:-1], analysis[N][1:]):
                try:
                    # FIXME: Better rate estimation formula!
                    # Hack to avoid division by zero
                    e0 = sqrt(a0["u_diff_norm"] / (a0["u_fine_norm"]+eps))
                    e1 = sqrt(a1["u_diff_norm"] / (a1["u_fine_norm"]+eps))
                    rate = e0 / (e1+eps)
                except:
                    rate = None
                rates.append(rate)
            print "Convergence rates w.r.t dt  at N=%g: %s" % (N, rates)


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

    return make_suite(TestGridConvergence, [schemes, problems])
