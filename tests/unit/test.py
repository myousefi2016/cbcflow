#!/usr/bin/env python
import io, os, sys, unittest, glob

def inject_python_path(modulepath):
    "Insert path into global python path, intended to import local module instead of installed one."
    if modulepath:
        sys.path.insert(0, modulepath)

def import_module_under_test(modulename):
    "Import module under test and print some info about it."
    mod = __import__(modulename)
    print("Running tests with %s version %s, date %s, imported from\n%s" % (
            modulename, mod.__version__, mod.__date__, mod.__file__))
    return mod

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

def run(modulename, modulepath, verbosity):
    "Import and run tests."
    # Import module
    inject_python_path(modulepath)
    mod = import_module_under_test(modulename)

    # Get tests
    tests = find_tests_in_current_dir()
    fullsuite = load_test_suite(tests)

    # Execute unittest framework
    runner = unittest.TextTestRunner(verbosity=verbosity)
    runner.run(fullsuite)

if __name__ == "__main__":
    args = sys.argv[1:]
    # TODO: Figure out a clean way to run both single and all tests locally or globally
    run(modulename="headflow",
        modulepath="../../site-packages",
        verbosity=2)

