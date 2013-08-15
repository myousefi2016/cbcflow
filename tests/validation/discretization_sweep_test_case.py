#!/usr/bin/env python

import sys, os, itertools, inspect
import unittest

from headflow import *

from math import sqrt
import dolfin
dolfin.parameters["allow_extrapolation"] = True


def run_discretization_sweep(scheme_factory, problem_factory, fields_factory, casedirs_factory, keys):
    "Call _run for each combination of N and dt and return dict with all data."
    data = {}
    for key in keys:
        # Make sure key is a tuple for convenience below
        if isinstance(key, list):
            key = tuple(key)
        elif not isinstance(key, tuple):
            key = (key,)

        # Construct fresh scheme
        scheme = scheme_factory(*key)

        # Construct fresh problem
        problem = problem_factory(*key)

        # Construct fresh fields
        fields = fields_factory()

        # Construct fresh casedirs
        casedirs = casedirs_factory(*((problem, scheme)+key))

        # Execute problem and collect postprocessing results
        results = run_problem(problem, scheme, fields, casedirs)

        # Store by parameter key
        data[key] = results

    return data

def run_problem(problem, scheme, fields, casedirs):
    "Execute problem with this scheme and return results collected from fields."

    # Construct fresh postprocessor
    casedir = os.path.join(*casedirs)
    pp = NSPostProcessor(ParamDict(casedir=casedir))

    # Add freshly constructed fields
    pp.add_fields(fields)

    # Construct the solver
    solver = NSSolver(problem, scheme, pp)

    # We allow just the solve to fail gracefully;
    # if other parts fail they can crash and burn
    try:
        # Execute solver (this may take some time...)
        ns = solver.solve()

        # Extract results for each field
        results = { f.name: pp.get(f.name) for f in fields }

        # Also store solver namespace
        results["namespace"] = ns

    except Exception as e:
        # No results, something failed so we just store the exception for the record
        results = { "exception": str(e) }

    return results


class DiscretizationSweepTestCase(unittest.TestCase):
    def init(self, scheme_factory, problem_factory, config={}):
        "Initialize this test case with given parameters."
        # These are callables to produce fresh NSScheme/NSProblem instances
        self._scheme_factory = scheme_factory
        self._problem_factory = problem_factory
        self.__Ns = config.get("Ns") if config else None
        self.__dts = config.get("dts") if config else None

    def shortDescription(self):
        "Returns a description of this test case for better reporting of failures in the unittest framework."
        sc = self._scheme_factory.func_code
        pc = self._problem_factory.func_code
        assert sc.co_filename == pc.co_filename
        loc = "%s:%d,%d" % (sc.co_filename, sc.co_firstlineno, pc.co_firstlineno)
        scode = inspect.getsource(self._scheme_factory)
        pcode = inspect.getsource(self._problem_factory)
        strings = (self.__class__.__name__, scode.lstrip(), pcode.lstrip(), loc)
        return "%s with params:\n    sf = %s    pf = %s    @ %s" % strings

    def _Ns(self):
        "Return range of spatial discretization parameters."
        if self.__Ns:
            return self.__Ns
        raise NotImplementedError("Missing implementation of _Ns in test class!")

    def _dts(self):
        "Return range of temporal discretization parameters."
        if self.__dts:
            return self.__dts
        raise NotImplementedError("Missing implementation of _dts in test class!")

    def _make_fields(self):
        "Return postprocessing fields to apply in solve."
        raise NotImplementedError("Missing implementation of _make_fields in test class!")

    def _analyse_data(self, data):
        "Analyse the data provided by the discretization parameter sweep."
        raise NotImplementedError("Missing implementation of _analyse_data in test class!")

    def _write_reference(self, name, data):
        path = os.path.join("output", self.__class__.__name__)
        os.makedirs(path)
        with open(os.path.join(path, name), "w") as f:
            f.write(data)

    def _read_reference(self, name):
        data = None
        path = os.path.join("references", self.__class__.__name__)
        with open(os.path.join(path, name), "r") as f:
            data = f.read()
        return data

    def runTest(self):
        "Run sweep over discretization parameters and delegate analysis to subclass"

        keys = [(N, dt) for N in self._Ns() for dt in self._dts()]

        scheme_factory = self._scheme_factory
        problem_factory = self._problem_factory
        fields_factory = self._make_fields

        casedirs_factory = lambda p, s, N, dt: [
            self.__class__.__name__,
            str(s),
            str(p),
            "N_%s" % str(N),
            "dt_%g" % dt,
            ]

        # Run simulations and collect all data
        data = run_discretization_sweep(scheme_factory, problem_factory, fields_factory, casedirs_factory, keys)

        # Analyse this data (implemented by subclass)
        self._analyse_data(data)


def make_suite(testclass, initargs):
    tests = []
    for args in itertools.product(*initargs):
        tc = testclass()
        tc.init(*args)
        tests.append(tc)
    ts = unittest.TestSuite(tests)
    return ts


def load_tests(*args):
    # Don't run the base class here as a test in auto-discover mode
    return []
