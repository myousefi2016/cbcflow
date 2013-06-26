#!/usr/bin/env python

import itertools
import unittest

def make_suite(testclass, initargs): # TODO: Move to a shared test utilities file
    tests = []
    for args in itertools.product(*initargs):
        tc = testclass()
        tc.init(*args)
        tests.append(tc)
    ts = unittest.TestSuite(tests)
    return ts

class TestTemplate(unittest.TestCase):
    "Example showing how to set up test cases parameterized by any number of parmeters."
    def init(self, parameter1, parameter2):
        self.parameter1 = parameter1
        self.parameter2 = parameter2

    def runTest(self):
        result = "Result of running some test with parameters %s and %s" % (self.parameter1, self.parameter2)
        self.assertIn(self.parameter1, result)
        self.assertIn(self.parameter2, result)

def load_tests(loader, standard_tests, none):
    parameters1 = ["one", "two"]
    parameters2 = ["abc", "def"]
    return make_suite(TestTemplate, [parameters1, parameters2])
