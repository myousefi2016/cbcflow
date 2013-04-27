#!/usr/bin/env python
"""
Tests of an alternative work in progress postprocessing framework in headflow.
"""

import unittest
from init_test import init_test
init_test(__name__)

from headflow import ParamDict
from headflow.core.parameterized import Parameterized

class NSPostProcessor2(Parameterized):
    def __init__(self, params=None):
        Parameterized.__init__(self, params)

        self._fields = {}
        self._cache = {}
        self._cache[0] = {}

    def add_field(self, field):
        assert field is self._fields.get(field.name, field)
        self._fields[field.name] = field

    def add_fields(self, fields):
        for field in fields:
            self.add_field(field)

    def get(self, name, timestep=0):
        c = self._cache[timestep]
        if name in c:
            v = c[name]
        else:
            f = self._fields[name]
            v = f.compute(self)
            c[name] = v
        return v

    def _need_for(self, action, name, t, timestep):
        f = self._fields[name]
        # TODO: Match f.params[action].* and t,timestep to see if we should do this now
        doit = False
        return doit

    def _do(self, action, name, t, timestep, value):
        if action == "save":
            pass#self._save(name, t, timestep, value)
        elif action == "plot":
            pass#self._plot(name, t, timestep, value)
        else:
            error("Unknown action %s." % action)

    def update_all(self, u, p, t, timestep):
        self._sorted_fields = self._fields.keys() # TODO: Make a topological ordering

        for name in self._sorted_fields:
            for action in ("save", "plot"):
                if self._need_for(action, name, t, timestep):
                    value = self.get(name)
                    self._do(action, name, t, timestep, value)

class PPField2(Parameterized):
    def __init__(self, params=None):
        Parameterized.__init__(self, params)

    @property
    def name(self):
        return self.__class__.__name__

    def compute(self, pp):
        raise NotImplementedError("A PPField must implement the compute function!")

class Velocity(PPField2):
    def __init__(self, params=None):
        PPField2.__init__(self, params)
        self.touched = 0

    def compute(self, pp):
        self.touched += 1
        return "u"

class Pressure(PPField2):
    def __init__(self, params=None):
        PPField2.__init__(self, params)
        self.touched = 0

    def compute(self, pp):
        self.touched += 1
        return "p"

class VelocityGradient(PPField2):
    def __init__(self, params=None):
        PPField2.__init__(self, params)
        self.touched = 0

    def compute(self, pp):
        self.touched += 1
        u = pp.get("Velocity")
        return "grad(%s)" % u

class Strain(PPField2):
    def __init__(self, params=None):
        PPField2.__init__(self, params)
        self.touched = 0

    def compute(self, pp):
        self.touched += 1
        Du = pp.get("VelocityGradient")
        return "epsilon(%s)" % Du

class Stress(PPField2):
    def __init__(self, params=None):
        PPField2.__init__(self, params)
        self.touched = 0

    def compute(self, pp):
        self.touched += 1
        epsilon = pp.get("Strain")
        p = pp.get("Pressure")
        return "sigma(%s, %s)" % (epsilon, p)

class TestPostProcessing2(unittest.TestCase):
    def test_(self):
        pp = NSPostProcessor2()
        u = Velocity()
        p = Pressure()
        Du = VelocityGradient()
        epsilon = Strain()
        sigma = Stress()
        pp.add_fields([u, p, Du, epsilon, sigma])

        self.assertEqual(u.touched, 0)

        self.assertEqual(pp.get("Strain"), "epsilon(grad(u))")
        self.assertEqual(u.touched, 1)
        self.assertEqual(Du.touched, 1)
        self.assertEqual(epsilon.touched, 1)
        self.assertEqual(p.touched, 0)
        self.assertEqual(sigma.touched, 0)

        self.assertEqual(pp.get("Stress"), "sigma(epsilon(grad(u)), p)")
        self.assertEqual(u.touched, 1) # Not recomputed!
        self.assertEqual(Du.touched, 1) # Not recomputed!
        self.assertEqual(epsilon.touched, 1) # Not recomputed!
        self.assertEqual(p.touched, 1)
        self.assertEqual(sigma.touched, 1)

if __name__ == "__main__":
    unittest.main()

