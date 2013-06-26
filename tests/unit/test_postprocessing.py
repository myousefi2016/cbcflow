#!/usr/bin/env python
"""
Tests of an alternative work in progress postprocessing framework in headflow.
"""

import unittest

from headflow import ParamDict, NSPostProcessor, WSS, Stress, Strain, Velocity, Pressure, VelocityGradient
from headflow.core.parameterized import Parameterized
'''
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
'''
'''
class PPField2(Parameterized):
    def __init__(self, params=None):
        Parameterized.__init__(self, params)

    @classmethod
    def default_base_params(cls):
        params = ParamDict(
            start_timestep = -1e16,
            end_timestep = 1e16,
            stride_timestep = 1,
            )
        return params

    @property
    def name(self):
        return self.__class__.__name__

    def init(self, pp): # TODO: Arguments?
        "Called prior to the simulation timeloop."
        pass

    def finalize(self, pp): # TODO: Arguments?
        "Called after the simulation timeloop."
        pass

    def compute(self, pp): # TODO: Arguments?
        "Called each time the quantity should be computed."
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

class TimeDerivative(PPField2):
    def __init__(self, name, params=None):
        PPField2.__init__(self, params)
        self.name = name

    def compute(self, pp):
        u1 = pp.get(self.name)
        u0 = pp.get(self.name, -1)

        t1 = pp.get("t")
        t0 = pp.get("t", -1)
        dt = dt-dt

        return (u1-u0)/dt
'''

class TestPostProcessing2(unittest.TestCase):
    def test_(self):
        pp = NSPostProcessor()
        u = Velocity()
        p = Pressure()
        Du = VelocityGradient()
        epsilon = Strain()
        sigma = Stress()
        pp.add_fields([u, p, Du, epsilon, sigma])
        
        # TODO: Figure out some proper tests for this.
        
        
        #self.assertEqual(u.touched, 0)

        #self.assertEqual(pp.get("Strain"), "epsilon(grad(u))")
        #self.assertEqual(u.touched, 1)
        #self.assertEqual(Du.touched, 1)
        #self.assertEqual(epsilon.touched, 1)
        #self.assertEqual(p.touched, 0)
        #self.assertEqual(sigma.touched, 0)

        #self.assertEqual(pp.get("Stress"), "sigma(epsilon(grad(u)), p)")
        #self.assertEqual(u.touched, 1) # Not recomputed!
        #self.assertEqual(Du.touched, 1) # Not recomputed!
        #self.assertEqual(epsilon.touched, 1) # Not recomputed!
        #self.assertEqual(p.touched, 1)
        #self.assertEqual(sigma.touched, 1)
        
