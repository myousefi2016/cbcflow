#!/usr/bin/env python
"""
Tests of postprocessing framework in headflow.
"""

import unittest

from collections import defaultdict

from headflow import ParamDict, NSProblem, NSPostProcessor, PPField, Velocity, Pressure, VelocityGradient, Strain, Stress, WSS
from headflow.core.spaces import NSSpacePoolSplit

import dolfin
from dolfin import UnitSquareMesh, Function, Expression, norm, errornorm, assemble, dx
from math import sqrt

# Avoid negative norms caused by instable tensor representation:
dolfin.parameters["form_compiler"]["representation"] = "quadrature"

_mesh = UnitSquareMesh(32, 32)

class MockProblem(NSProblem):
    def __init__(self):
        NSProblem.__init__(self)
        self.initialize_geometry(_mesh)

    @classmethod
    def default_user_params(cls):
        return ParamDict(
            mu = 0.1,
            rho = 0.9,
            )

class MockPPField(PPField):
    def __init__(self, params=None):
        PPField.__init__(self, params)
        self.touched = 0

    @classmethod
    def default_user_params(cls):
        return ParamDict(
            # Don't compute unless asked
            start_timestep=1e16,
            end_timestep=-1e16,
            stride_timestep=0,

            start_time=1.0e16,
            end_time=-1.0e16,
            stride_time=0.0,

            # Don't save or plot, but call callback
            save = False,
            plot = False,
            callback = True,
            )


class MockVelocity(MockPPField):
    def __init__(self, params=None):
        MockPPField.__init__(self, params)

    def compute(self, pp, spaces, problem):
        self.touched += 1
        return "u"

class MockPressure(MockPPField):
    def __init__(self, params=None):
        MockPPField.__init__(self, params)

    def compute(self, pp, spaces, problem):
        self.touched += 1
        return "p"

class MockVelocityGradient(MockPPField):
    def __init__(self, params=None):
        MockPPField.__init__(self, params)

    def compute(self, pp, spaces, problem):
        self.touched += 1
        u = pp.get("MockVelocity")
        return "grad(%s)" % u

class MockStrain(MockPPField):
    def __init__(self, params=None):
        MockPPField.__init__(self, params)

    def compute(self, pp, spaces, problem):
        self.touched += 1
        Du = pp.get("MockVelocityGradient")
        return "epsilon(%s)" % Du

class MockStress(MockPPField):
    def __init__(self, params=None):
        MockPPField.__init__(self, params)

    def compute(self, pp, spaces, problem):
        self.touched += 1
        epsilon = pp.get("MockStrain")
        p = pp.get("MockPressure")
        return "sigma(%s, %s)" % (epsilon, p)

class MockTimeDerivative(MockPPField):
    def __init__(self, name, params=None):
        MockPPField.__init__(self, params)
        self.name = name

    def compute(self, pp, spaces, problem):
        self.touched += 1

        u1 = pp.get(self.name)
        u0 = pp.get(self.name, -1)

        t1 = pp.get("t")
        t0 = pp.get("t", -1)
        dt = t1-t0

        return "(%s - %s) / (%g)" % (u1, u0, dt)


class TestPostProcessing2(unittest.TestCase):
    def test_mock_fields_get_correct_compute_calls_single_timestep(self):
        # Create mock field objects
        pressure = MockPressure()
        velocity = MockVelocity()
        Du = MockVelocityGradient()
        epsilon = MockStrain()
        sigma = MockStress()

        # Add fields to postprocessor
        pp = NSPostProcessor()
        pp.add_fields([pressure, velocity, Du, epsilon, sigma])
        # Adding in random ordering not respecting dependencies
        # results in an assertion failure in nspostprocessor
        # because the type of the built-in Velocity is not the
        # same as MockVelocity, this is ok I guess.
        #pp.add_fields([Du, pressure, velocity, sigma, Du, epsilon])

        # Nothing has been computed yet
        self.assertEqual(velocity.touched, 0)
        self.assertEqual(Du.touched, 0)
        self.assertEqual(epsilon.touched, 0)
        self.assertEqual(pressure.touched, 0)
        self.assertEqual(sigma.touched, 0)

        # Mock problem
        problem = MockProblem()

        # Mock internal scheme state
        spaces = NSSpacePoolSplit(problem.mesh, 1, 1)
        u = Function(spaces.V)
        u.interpolate(Expression(("x[0]", "2*x[1]")))
        p = Function(spaces.Q)
        p.interpolate(Expression("3*x[0]*x[1]"))
        t = 0.0
        timestep = 0

        # Update postprocessor fields using mock problem
        pp.update_all(u, p, t, timestep, spaces, problem)

        # Get strain twice
        for i in range(2):
            strain = pp.get("MockStrain")
            # Check value
            self.assertEqual(strain, "epsilon(grad(u))")
            # Check the right things are computed but only the first time
            self.assertEqual(velocity.touched, 1) # Only increased first iteration!
            self.assertEqual(Du.touched, 1) # ...
            self.assertEqual(epsilon.touched, 1) # ...
            self.assertEqual(pressure.touched, 0) # Not computed!
            self.assertEqual(sigma.touched, 0) # ...

        # Get stress twice
        for i in range(2):
            stress = pp.get("MockStress")
            # Check value
            self.assertEqual(stress, "sigma(epsilon(grad(u)), p)")
            # Check the right things are computed but only the first time
            self.assertEqual(velocity.touched, 1) # Not recomputed!
            self.assertEqual(Du.touched, 1) # ...
            self.assertEqual(epsilon.touched, 1) # ...
            self.assertEqual(pressure.touched, 1) # Only increased first iteration!
            self.assertEqual(sigma.touched, 1) # ...

    def test_get_time_and_solution(self):

        # This is the object we want to test!
        pp = NSPostProcessor()

        # Attach a callback to postprocessor so we can inspect direct compute requests
        def ppcallback(field, data, t, timestep):
            ppcallback.calls[field.name].append((t, timestep))
        ppcallback.calls = defaultdict(list)
        pp._callback = ppcallback

        # Setup some mock problem state
        problem = MockProblem()

        # Toy state evolution
        uexpr = Expression(("2.0 + t + x[0]", "3.0 + t * 2 * x[1]"), t=0.0)
        pexpr = Expression("1.0 + t * 3 * x[0] * x[1]", t=0.0)

        # Setup some mock scheme state
        spaces = NSSpacePoolSplit(problem.mesh, 1, 1)
        u = Function(spaces.V)
        p = Function(spaces.Q)

        # Set initial time
        t = 0.0
        timestep = 0
        dt = 0.1

        # Set "initial condition"
        uexpr.t = t
        pexpr.t = t
        u.interpolate(uexpr)
        p.interpolate(pexpr)

        # Update postprocessor, this is where the main code under test is
        pp.update_all(u, p, t, timestep, spaces, problem)

        # Check that we recover the velocity, pressure, and time
        self.assertEqual(t, pp.get("t"))
        self.assertEqual(timestep, pp.get("timestep"))
        self.assertAlmostEqual(norm(u), norm(pp.get("Velocity")))
        self.assertAlmostEqual(norm(p), norm(pp.get("Pressure")))
        self.assertGreater(norm(u), 0.0) # Check that previous test was not worthless
        self.assertGreater(norm(p), 0.0)

        # Loop over "timesteps"
        for k in range(4):
            t += dt
            timestep += 1

            # Update "solution"
            uexpr.t = t
            pexpr.t = t
            u.interpolate(uexpr)
            p.interpolate(pexpr)

            # Update postprocessor, this is where the main code under test is
            pp.update_all(u, p, t, timestep, spaces, problem)

            # Check that we recover the basic quantities velocity, pressure, and time
            self.assertEqual(t, pp.get("t"))
            self.assertEqual(timestep, pp.get("timestep"))
            self.assertAlmostEqual(norm(u), norm(pp.get("Velocity")))
            self.assertAlmostEqual(norm(p), norm(pp.get("Pressure")))
            self.assertGreater(norm(u), 0.0) # Check that previous test was not worthless
            self.assertGreater(norm(p), 0.0)

        # We didn't make any direct compute requests, so no actions should be triggered:
        self.assertEqual(ppcallback.calls, defaultdict(list))

    def test_get_compute_velocity_gradient_strain_and_stress(self):
        # FIXME: Test these computations with both mixed, split and segregated "schemes"!

        # This is the object we want to test!
        pp = NSPostProcessor()
        pp.add_fields([VelocityGradient(), Pressure(), Stress(), Velocity(), Strain()])

        # Attach a callback to postprocessor so we can inspect direct compute requests
        def ppcallback(field, data, t, timestep):
            ppcallback.calls[field.name].append((t, timestep))
        ppcallback.calls = defaultdict(list)
        pp._callback = ppcallback

        # Setup some mock problem state
        problem = MockProblem()

        # Toy state evolution
        uexpr = Expression(("1.0 + 2*x[0] + 3*x[1]", "1.0 + 4*x[0] + 5*x[1]"))
        pexpr = Expression("1.0 + 3 * x[0] * x[1]")

        # Manually derived quantities
        Du_expr = Expression((("2.0", "3.0"), ("4.0", "5.0")))
        epsilon_expr = Expression((("2.0", "3.5"), ("3.5", "5.0")))
        sigma_expr = Expression((("2*0.2 - (1.0 + 3 * x[0] * x[1])", "2*0.35"),
                             ("2*0.35", "2*0.5 - (1.0 + 3 * x[0] * x[1])")))
        self.assertAlmostEqual(problem.params.mu, 0.1) # Using this in sigma_expr

        # Setup some mock scheme state
        spaces = NSSpacePoolSplit(problem.mesh, 1, 1)
        u = Function(spaces.V)
        p = Function(spaces.Q)

        # Set initial time
        t = 0.0
        timestep = 0
        dt = 0.1

        # Set "initial condition"
        u.interpolate(uexpr)
        p.interpolate(pexpr)

        # Update postprocessor, this is where the main code under test is
        pp.update_all(u, p, t, timestep, spaces, problem)

        # Check that we recover the velocity, pressure, and time
        def _errornorm(expr, name):
            a = assemble(expr**2*dx(), mesh=problem.mesh)
            e = assemble(pp.get(name)**2*dx(), mesh=problem.mesh)
            d = assemble((expr - pp.get(name))**2*dx(), mesh=problem.mesh)
            print a, e, d, d/e
            return d / e

        if 0:
            print "u", _errornorm(uexpr, "Velocity")
            print "p", _errornorm(pexpr, "Pressure")
            print "g", _errornorm(Du_expr, "VelocityGradient")
            print "s", _errornorm(epsilon_expr, "Strain")

        self.assertAlmostEqual(_errornorm(uexpr, "Velocity"), 0.0)
        self.assertAlmostEqual(_errornorm(pexpr, "Pressure"), 0.0)

        if 0:
            dolfin.plot(Du_expr[0,:], mesh=problem.mesh)
            dolfin.plot(Du_expr[1,:], mesh=problem.mesh)
            dolfin.plot(pp.get("VelocityGradient")[0,:], mesh=problem.mesh)
            dolfin.plot(pp.get("VelocityGradient")[0,:], mesh=problem.mesh)
            dolfin.interactive()

        self.assertAlmostEqual(_errornorm(Du_expr, "VelocityGradient"), 0.0)
        self.assertAlmostEqual(_errornorm(epsilon_expr, "Strain"), 0.0)
        self.assertLess(_errornorm(sigma_expr, "Stress"), 3e-4) # This is less accurate, but converges

        # We didn't make any direct compute requests, so no actions should be triggered:
        self.assertEqual(ppcallback.calls, defaultdict(list))


    def test_get_wss(self):
        self.assertEqual(1, 1) # FIXME

    def test_get_first_time_derivative(self):
        self.assertEqual(1, 1) # FIXME

    def test_get_second_time_derivative(self):
        self.assertEqual(1, 1) # FIXME

    def test_get_third_time_derivative(self):
        self.assertEqual(1, 1) # FIXME
