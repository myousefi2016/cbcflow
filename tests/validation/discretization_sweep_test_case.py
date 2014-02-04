#!/usr/bin/env python

import sys, os, itertools, inspect
import unittest

from cbcflow import *

from math import sqrt
import dolfin
dolfin.parameters["allow_extrapolation"] = True


def run_problem(problem, scheme, fields, casedir):
    "Execute problem with this scheme and return results collected from fields."

    # Construct fresh postprocessor
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

    except Exception:
        # No results, something failed so we just store the exception for the record
        import sys, traceback
        ex_type, ex, tb = sys.exc_info()
        tb = ''.join(traceback.format_tb(tb))
        results = { "exception": ex, "traceback": tb }

    return results


def run_discretization_sweep(scheme_factory, problem_factory, fields_factory, casedirs_factory, keys):
    """Call _run for each combination of N and dt and return dict with all data.

    TODO: Document factory methods.
    """
    all_casedirs = set()
    data = {}
    for key in keys:
        # Make sure key is a tuple for convenience below
        if isinstance(key, list):
            key = tuple(key)
        elif not isinstance(key, tuple):
            key = (key,)

        # Construct fresh scheme
        scheme = scheme_factory()

        # Construct fresh problem
        problem = problem_factory(*key)

        # Construct fresh fields
        fields = fields_factory()

        # Construct fresh casedirs
        casedirs = tuple(casedirs_factory(*((problem, scheme)+key)))
        all_casedirs.add(casedirs)
        casedir = os.path.join(*casedirs)

        # Execute problem and collect postprocessing results
        results = run_problem(problem, scheme, fields, casedir)

        # Store by parameter key
        data[key] = results

    # Extract common part of casedirs
    data["basedir"] = common_basedir(all_casedirs)

    return data


def common_basedir(all_casedirs):
    all_casedirs = list(all_casedirs)
    n = len(all_casedirs)
    j = min(len(a) for a in all_casedirs)

    for k0 in xrange(n):
        a = all_casedirs[k0]
        for k1 in xrange(k0+1, n):
            b = all_casedirs[k1]

            for i in xrange(min(j,len(b))):
                if a[i] != b[i]:
                    j = min(j, i)
                    break
    basedir = all_casedirs[0][:j]
    assert all(basedir == cn[:j] for cn in all_casedirs)
    return basedir


class DiscretizationSweepTestCase(unittest.TestCase):
    def init(self, scheme_factory, problem_factory, config={}):
        "Initialize this test case with given parameters."
        # These are callables to produce fresh NSScheme/NSProblem instances
        self._scheme_factory = scheme_factory
        self._problem_factory = problem_factory
        self.__refinement_levels = config.get("refinement_levels") if config else None
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

    def _refinement_levels(self):
        "Return range of spatial discretization parameters."
        if self.__refinement_levels:
            return self.__refinement_levels
        raise NotImplementedError("Missing implementation of _refinement_levels in test class!")

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

    def _write_reference(self, path, name, data, basedir="output"):
        path = (path,) if isinstance(path, str) else tuple(path)
        path = os.path.join(*((basedir,) + path))
        print "PATH", path
        if not os.path.isdir(path):
            os.makedirs(path)
        fullname = os.path.join(path, name)
        with open(fullname, "w") as f:
            f.write(data)

    def _read_reference(self, path, name, basedir="references"):
        data = None
        path = (path,) if isinstance(path, str) else tuple(path)
        path = os.path.join(*((basedir,) + path))
        fullname = os.path.join(path, name)
        if os.path.isfile(fullname):
            with open(fullname, "r") as f:
                data = f.read()
        else:
            data = None
        return data

    def runTest(self):
        "Run sweep over discretization parameters and delegate analysis to subclass"

        keys = [(refinement_level, dt) for refinement_level in self._refinement_levels() for dt in self._dts()]

        scheme_factory = self._scheme_factory
        problem_factory = self._problem_factory
        fields_factory = self._make_fields

        casedirs_factory = lambda p, s, refinement_level, dt: [
            self.__class__.__name__,
            s.shortname(),
            p.shortname(),
            "ref_%s" % str(refinement_level),
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
