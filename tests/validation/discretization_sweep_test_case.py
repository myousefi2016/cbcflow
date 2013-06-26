#!/usr/bin/env python

import sys, os, itertools, inspect
import unittest

from headflow import *

from math import sqrt
import dolfin
dolfin.parameters["allow_extrapolation"] = True


class DiscretizationSweepTestCase(unittest.TestCase):
    def init(self, scheme_factory, problem_factory, config={}):
        "Initialize this test case with given parameters."
        # These are callables to produce fresh NSScheme/NSProblem instances
        self.sf = scheme_factory
        self.pf = problem_factory
        self.__Ns = config.get("Ns") if config else None
        self.__dts = config.get("dts") if config else None

    def shortDescription(self):
        sc = self.sf.func_code
        pc = self.pf.func_code
        assert sc.co_filename == pc.co_filename
        loc = "%s:%d,%d" % (sc.co_filename, sc.co_firstlineno, pc.co_firstlineno)
        scode = inspect.getsource(self.sf)
        pcode = inspect.getsource(self.pf)

        return "%s with params:\n    sf = %s    pf = %s    @ %s" % (self.__class__.__name__, scode.lstrip(), pcode.lstrip(), loc)

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

    def _run_discretization_sweep(self, Ns, dts):
        "Call _run for each combination of N and dt and return dict with all data."
        data = {}
        for N in Ns:
            for dt in dts:
                data[(N,dt)] = self._run(N, dt)
        return data

    def _run(self, N, dt):
        "Execute problem for this N,dt and collect results from fields."

        # Construct fresh scheme
        s = self.sf()

        # Construct fresh problem
        p = self.pf(N, dt)

        # Construct fresh postprocessor
        dirs = [self.__class__.__name__, str(s), str(p), "N_%s"%str(N), "dt_%g"%dt]
        casedir = os.path.join(*dirs)
        pp = NSPostProcessor(ParamDict(casedir=casedir))

        # Add freshly constructed fields
        fields = self._make_fields()
        pp.add_fields(fields)

        # Construct the solver
        solver = NSSolver(p, s, pp)

        # We allow just the solve to fail gracefully;
        # if other parts fail they can crash and burn
        try:
            ns = solver.solve()

            # Extract results
            # TODO: Standardize and document interface for ppfield.get_data()
            #results = { f.__class__.__name__: f.get_data() for f in fields }
            results = {f.name: pp._cache[0][f.name] for f in fields}
            results["namespace"] = ns
        except Exception as e:
            # No results
            
            results = { "exception": str(e) }

        return results

    def runTest(self):
        "Run sweep over discretization parameters and delegate analysis to subclass"
        # Run simulations and collect all data (implemented above)
        data = self._run_discretization_sweep(self._Ns(), self._dts())
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
