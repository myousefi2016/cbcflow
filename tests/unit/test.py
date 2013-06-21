#!/usr/bin/env python

from init_test import init_test
init_test(__name__)

import os, sys, unittest, glob

def inject_python_path(modulepath):
    "Insert path into global python path, intended to import local module instead of installed one."
    if modulepath:
        sys.path.insert(0, modulepath)

def find_tests_in_current_dir():
    "Find test_*.py in current directory and return basenames."
    testpattern = "test_*.py"
    tests = sorted(os.path.basename(f).replace(".py", "") for f in glob.glob(testpattern))
    return tests

def load_test_suite(tests):
    "Build unittest suite from list of importable module names."
    loader = unittest.TestLoader()
    fullsuite = unittest.TestSuite()
    for case in tests:
        casemodule = __import__(case)
        casesuite = loader.loadTestsFromModule(casemodule)
        fullsuite.addTests(casesuite)
    return fullsuite

def discover_and_run_tests(modulename, verbosity):
    "Import and run tests."
    # Get tests
    tests = find_tests_in_current_dir()
    fullsuite = load_test_suite(tests)

    # Execute unittest framework
    runner = unittest.TextTestRunner(verbosity=verbosity)
    results = runner.run(fullsuite)
    return len(results.errors)+len(results.failures)

if __name__ == "__main__":
    num_fails = discover_and_run_tests(modulename="headflow", verbosity=2)
    sys.exit(num_fails)
