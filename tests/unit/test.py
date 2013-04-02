#!/usr/bin/env python
import io, os, sys, unittest, glob

def run(modulename, modulepath, verbosity):
    if modulepath is not None:
        # Insert path, intended to import local module instead of installed one
        sys.path.insert(0, modulepath)
    mod = __import__(modulename)
    print("Running tests with %s version %s, date %s." % (
            modulename, mod.__version__, mod.__date__))

    # Running tests from all test_foo.py files
    testpattern = "test_*.py"
    tests = sorted(os.path.basename(f).replace(".py", "") for f in glob.glob(testpattern))

    # Setup unittest 
    loader = unittest.TestLoader()
    fullsuite = unittest.TestSuite()
    for case in tests:
        casemodule = __import__(case)
        casesuite = loader.loadTestsFromModule(casemodule)
        fullsuite.addTests(casesuite)
    unittest.TextTestRunner(verbosity=verbosity).run(fullsuite)

if __name__ == "__main__":
    run(modulename="headflow",
        modulepath="../../site-packages",
        verbosity=2)
