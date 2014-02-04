#!/usr/bin/env python

import sys, os, itertools
import unittest

from cbcflow import *

from math import sqrt
import dolfin
from dolfin import *
dolfin.parameters["allow_extrapolation"] = True

from discretization_sweep_test_case import DiscretizationSweepTestCase, make_suite

_logfile = open("testlog.log","w")

class TestGridConvergence(DiscretizationSweepTestCase):

    def _refinement_levels(self):
        "Return range of spatial discretization parameters."
        return [0,1]

    def _dts(self):
        "Return range of temporal discretization parameters."
        return [0.05, 0.025]

    def _make_fields(self):
        "Return postprocessing fields to apply in solve."
        return [Velocity()] # FIXME: Use VelocityError(), PressureError() here?

    def _analyse_data(self, data):
        "Analyse the data provided by the discretization parameter sweep."
        self._analyze_spatial_convergence(data)
        self._analyze_temporal_convergence(data)

    def _print(self, msg):
        # FIXME: Don't use print! Log to file and compare with reference data,
        #        and use assertions for values with known properties.
        print msg
        _logfile.write(msg + "\n")

    def _analyze_spatial_convergence(self, data):
        """With time discretization fixed, loop over
        spatial discretization and compute norms."""

        self._print("BEGIN _analyze_spatial_convergence")

        eps = 1e-14

        refinement_levels = self._refinement_levels()
        dts = self._dts()

        analysis = {}
        for dt in dts:
            analysis[dt] = []
            for refinement_level0, refinement_level1 in zip(refinement_levels[:-1], refinement_levels[1:]):
                # Fetch data at two refinements
                key0 = (refinement_level0, dt)
                key1 = (refinement_level1, dt)
                data0 = data[key0]
                data1 = data[key1]

                # Validate data
                ex = data0.get("exception")
                if ex is not None:
                    key = key0
                    tb = data0.get("traceback")
                    msg = "At refinement_level={refinement_level}; dt={dt}; got exception {ext}:\n{exs}\nTraceback:\n{tb}".format(refinement_level=key[0], dt=key[1], ext=type(ex), exs=str(ex), tb=tb)
                    self._print(msg)
                    continue
                ex = data1.get("exception")
                if ex is not None:
                    key = key1
                    tb = data1.get("traceback")
                    msg = "At refinement_level={refinement_level}; dt={dt}; got exception {ext}:\n{exs}\nTraceback:\n{tb}".format(refinement_level=key[0], dt=key[1], ext=type(ex), exs=str(ex), tb=tb)
                    self._print(msg)
                    continue

                # Do some analysis
                u_coarse = data0["Velocity"]
                u_fine = data1["Velocity"]
                uc = interpolate(u_coarse, u_fine.function_space())
                u_coarse_norm = assemble(u_coarse**2*dx())
                u_fine_norm = assemble(u_fine**2*dx())
                u_diff_norm = assemble((u_fine - uc)**2*dx())
                u_rel_norm = sqrt(u_diff_norm / (u_fine_norm+eps))

                # Package nicely
                a = {
                    "refinement_level0": key0[0],
                    "refinement_level1": key1[1],
                    "dt0": key0[1],
                    "dt1": key1[1],
                    "u_coarse_norm": u_coarse_norm,
                    "u_fine_norm": u_fine_norm,
                    "u_diff_norm": u_diff_norm,
                    "u_rel_norm": u_rel_norm,
                    }
                analysis[dt].append(a)

        # Slight modification of the prints that Kent did before:
        for key in sorted(analysis.keys()):
            for a in analysis[key]:
                msg = "At refinement_level0,dt0={refinement_level0},{dt0}; refinement_level1,dt1={refinement_level1},{dt1}; abs diff = {udn}; fine norm = {ufn}; relative norm = {urn}.".format(refinement_level0=a["refinement_level0"], refinement_level1=a["refinement_level1"],
                                                                                                                                    dt0=a["dt0"], dt1=a["dt1"],
                                                                                                                                    udn=a["u_diff_norm"],
                                                                                                                                    ufn=a["u_fine_norm"],
                                                                                                                                    urn=a["u_rel_norm"])
                self._print(msg)

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
            msg = "Convergence rates w.r.t N at dt=%g: %s" % (dt, rates)
            self._print(msg)

        self._print("END _analyze_spatial_convergence")

    def _analyze_temporal_convergence(self, data):
        "With spatial discretization fixed, loop over time discretization and compute norms."
        self._print("BEGIN _analyze_temporal_convergence")

        eps = 1e-14

        refinement_levels = self._refinement_levels()
        dts = self._dts()

        analysis = {}
        for refinement_level in refinement_levels:
            analysis[refinement_level] = []
            for dt0, dt1 in zip(dts[:-1], dts[1:]):
                # Fetch data at two refinements
                key0 = (refinement_level, dt0)
                key1 = (refinement_level, dt1)
                data0 = data[key0]
                data1 = data[key1]

                # Validate data
                ex = data0.get("exception")
                if ex is not None:
                    key = key0
                    tb = data0.get("traceback")
                    msg = "At refinement_level={refinement_level}; dt={dt}; got exception {ext}:\n{exs}\nTraceback:\n{tb}".format(refinement_level=key[0], dt=key[1], ext=type(ex), exs=str(ex), tb=tb)
                    self._print(msg)
                    continue
                ex = data1.get("exception")
                if ex is not None:
                    key = key1
                    tb = data1.get("traceback")
                    msg = "At refinement_level={refinement_level}; dt={dt}; got exception {ext}:\n{exs}\nTraceback:\n{tb}".format(refinement_level=key[0], dt=key[1], ext=type(ex), exs=str(ex), tb=tb)
                    self._print(msg)
                    continue

                # Do some analysis
                u_coarse = data0["Velocity"]
                u_fine = data1["Velocity"]
                uc = dolfin.interpolate(u_coarse, u_fine.function_space())
                u_coarse_norm = assemble(u_coarse**2*dx())
                u_fine_norm = assemble(u_fine**2*dx())
                u_diff_norm = assemble((u_fine - uc)**2*dx())
                u_rel_norm = sqrt(u_diff_norm / (u_fine_norm+eps))

                # Package nicely
                a = {
                    "refinement_level0": key0[0],
                    "refinement_level1": key1[1],
                    "dt0": key0[1],
                    "dt1": key1[1],
                    "u_coarse_norm": u_coarse_norm,
                    "u_fine_norm": u_fine_norm,
                    "u_diff_norm": u_diff_norm,
                    "u_rel_norm": u_rel_norm,
                    }
                analysis[refinement_level].append(a)

        # Slight modification of the prints that Kent did before:
        for key in sorted(analysis.keys()):
            for a in analysis[key]:
                msg = "At refinement_level0,dt0={refinement_level0},{dt0}; refinement_level1,dt1={refinement_level1},{dt1}; abs diff = {udn}; fine norm = {ufn}; relative norm = {urn}.".format(refinement_level0=a["refinement_level0"], refinement_level1=a["refinement_level1"],
                                                                                                                                    dt0=a["dt0"], dt1=a["dt1"],
                                                                                                                                    udn=a["u_diff_norm"],
                                                                                                                                    ufn=a["u_fine_norm"],
                                                                                                                                    urn=a["u_rel_norm"])
                self._print(msg)

        # FIXME: Add assertions to validate convergence rate from analysis[N][:]
        for refinement_level in refinement_levels:
            rates = []
            for a0, a1 in zip(analysis[refinement_level][:-1], analysis[refinement_level][1:]):
                try:
                    # FIXME: Better rate estimation formula!
                    # Hack to avoid division by zero
                    e0 = sqrt(a0["u_diff_norm"] / (a0["u_fine_norm"]+eps))
                    e1 = sqrt(a1["u_diff_norm"] / (a1["u_fine_norm"]+eps))
                    rate = e0 / (e1+eps)
                except:
                    rate = None
                rates.append(rate)
            msg = "Convergence rates w.r.t dt  at refinement_level=%g: %s" % (refinement_level, rates)
            self._print(msg)

        self._print("END _analyze_temporal_convergence")


# FIXME: Need a better solution for importing demos!
# Importing problems from cbcflow/demo/
# NB! Assuming run from the cbcflow/tests/ directory!
sys.path.insert(0, "../demo/undocumented/FlowAroundCylinder")
sys.path.insert(0, "../demo/undocumented/Beltrami")
from FlowAroundCylinder import FlowAroundCylinder
from Beltrami import Beltrami


def load_tests(loader, standard_tests, none):

    # FIXME: Make fast and slow suite

    # FIXME: Add more schemes, use all of official_schemes
    schemes = [
        lambda: IPCS(),
        lambda: IPCS_Stable(),
        lambda: IPCS_Stabilized({'theta':0.5}),
        ]

    # FIXME: Add more problems
    problems = [
        lambda refinement_level,dt: Beltrami(ParamDict(refinement_level=refinement_level, dt=dt, T=dt*2)), # FIXME: Limiting T for debugging
        lambda refinement_level,dt: FlowAroundCylinder(ParamDict(refinement_level=refinement_level, dt=dt, T=dt*2)), # FIXME: Limiting T for debugging
        ]

    return make_suite(TestGridConvergence, [schemes, problems])
