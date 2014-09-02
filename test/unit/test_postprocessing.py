#!/usr/bin/env py.test
"""
Tests of postprocessing framework in cbcflow.
"""

from collections import defaultdict

from cbcflow import (ParamDict, NSProblem, NSPostProcessor,
    PPField, Velocity, Pressure, VelocityGradient, Strain, Stress, WSS,
    TimeDerivative, SecondTimeDerivative, TimeIntegral, Norm)
from cbcflow.utils.core import NSSpacePoolSplit

import dolfin
from dolfin import UnitSquareMesh, Function, Expression, norm, errornorm, assemble, dx
from math import sqrt

# Avoid negative norms caused by instable tensor representation:
dolfin.parameters["form_compiler"]["representation"] = "quadrature"

class MockProblem(NSProblem):
    def __init__(self, params=None):
        NSProblem.__init__(self, params)
        mesh = UnitSquareMesh(32, 32)
        self.initialize_geometry(mesh)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            T = 1.0,
            dt = 0.1,
            mu = 0.1,
            rho = 0.9,
            )
        return params

ppf_delayed_cb_params = ParamDict(
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

ppf_immediate_cb_params = ParamDict(
    # Compute every timestep
    start_timestep=-1e16,
    end_timestep=1e16,
    stride_timestep=1,

    start_time=-1.0e16,
    end_time=1.0e16,
    stride_time=0.0,

    # Don't save or plot, but call callback
    save = False,
    plot = False,
    callback = True,
    
    # Return on compute
    finalize = False,
    )

class MockPPField(PPField):
    def __init__(self, params=None):
        PPField.__init__(self, params)
        self.touched = 0

    @classmethod
    def default_params(cls):
        params = PPField.default_params()
        params.replace(**ppf_delayed_cb_params)
        return params


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


def test_mock_fields_get_correct_compute_calls_single_timestep():
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
    assert velocity.touched == 0
    assert Du.touched == 0
    assert epsilon.touched == 0
    assert pressure.touched == 0
    assert sigma.touched == 0

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
    pp.update_all({"Velocity": u, "Pressure": p}, t, timestep, spaces, problem)

    # Get strain twice
    for i in range(2):
        strain = pp.get("MockStrain")
        # Check value
        assert strain == "epsilon(grad(u))"
        # Check the right things are computed but only the first time
        assert velocity.touched == 1 # Only increased first iteration!
        assert Du.touched == 1 # ...
        assert epsilon.touched == 1 # ...
        assert pressure.touched == 0 # Not computed!
        assert sigma.touched == 0 # ...

    # Get stress twice
    for i in range(2):
        stress = pp.get("MockStress")
        # Check value
        assert stress == "sigma(epsilon(grad(u)), p)"
        # Check the right things are computed but only the first time
        assert velocity.touched == 1 # Not recomputed!
        assert Du.touched == 1 # ...
        assert epsilon.touched == 1 # ...
        assert pressure.touched == 1 # Only increased first iteration!
        assert sigma.touched == 1 # ...

def test_get_time_and_solution_from_pp_after_update_all():

    # This is the object we want to test!
    pp = NSPostProcessor()

    # Should add Velocity and Pressure manually to make these available,
    # this is to be able to skip unnecessary copying of these functions in the pp
    ppp = MockPPField.default_params()
    pp.add_fields([Velocity(ppp), Pressure(ppp)])

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
    pp.update_all({"Velocity": u, "Pressure": p}, t, timestep, spaces, problem)

    # Check that we recover the time and timestep, these are not fields but handled specially
    assert t == pp.get("t")
    assert timestep == pp.get("timestep")
    # Check that we recover the velocity and pressure, these were requested above
    assert abs( (norm(u)) - (norm(pp.get("Velocity"))) ) < 1e-8
    assert abs( (norm(p)) - (norm(pp.get("Pressure"))) ) < 1e-8
    assert norm(u) > 0.0 # Check that previous test was not worthless
    assert norm(p) > 0.0

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
        pp.update_all({"Velocity": u, "Pressure": p}, t, timestep, spaces, problem)

        # Check that we recover the basic quantities velocity, pressure, and time
        assert t == pp.get("t")
        assert timestep == pp.get("timestep")
        assert abs( (norm(u)) - (norm(pp.get("Velocity"))) ) < 1e-8
        assert abs( (norm(p)) - (norm(pp.get("Pressure"))) ) < 1e-8
        assert norm(u) > 0.0 # Check that previous test was not worthless
        assert norm(p) > 0.0

    # We didn't make any direct compute requests, so no actions should be triggered:
    assert ppcallback.calls == defaultdict(list)

def test_get_compute_velocity_gradient_strain_and_stress():
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
    assert abs( (problem.params.mu) - (0.1) ) < 1e-8 # Using this in sigma_expr

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
    pp.update_all({"Velocity": u, "Pressure": p}, t, timestep, spaces, problem)

    # Check that we recover the velocity, pressure, and time
    def _errornorm(expr, name):
        a = assemble(expr**2*dx(), mesh=problem.mesh)
        e = assemble(pp.get(name)**2*dx(), mesh=problem.mesh)
        d = assemble((expr - pp.get(name))**2*dx(), mesh=problem.mesh)
        #print a, e, d, d/e
        return d / e

    if 0:
        print "u", _errornorm(uexpr, "Velocity")
        print "p", _errornorm(pexpr, "Pressure")
        print "g", _errornorm(Du_expr, "VelocityGradient")
        print "s", _errornorm(epsilon_expr, "Strain")

    assert _errornorm(uexpr, "Velocity") < 1e-8
    assert _errornorm(pexpr, "Pressure") < 1e-8

    if 0:
        dolfin.plot(Du_expr[0,:], mesh=problem.mesh)
        dolfin.plot(Du_expr[1,:], mesh=problem.mesh)
        dolfin.plot(pp.get("VelocityGradient")[0,:], mesh=problem.mesh)
        dolfin.plot(pp.get("VelocityGradient")[0,:], mesh=problem.mesh)
        dolfin.interactive()

    assert _errornorm(Du_expr, "VelocityGradient") < 1e-8
    assert _errornorm(epsilon_expr, "Strain") < 1e-8
    assert _errornorm(sigma_expr, "Stress") < 3e-4 # This is less accurate, but converges

    # We didn't make any direct compute requests, so no actions should be triggered:
    assert ppcallback.calls == defaultdict(list)

def test_get_wss():
    assert 1 == 1 # FIXME

def test_get_first_time_derivative():
    # This is the object we want to test!
    pp = NSPostProcessor()
    ppp = ppf_immediate_cb_params # Important that we ask for the fields to be computed, as time dependencies will not work retrospectively after time loop is over!
    pp.add_fields([
            #Velocity(ppp),
            Norm("Pressure", ppp),
            TimeDerivative("t", ppp),
            TimeDerivative("timestep", ppp),
            TimeDerivative("Norm_Pressure", ppp),
            TimeDerivative("Pressure", ppp),
            SecondTimeDerivative("t", ppp),
            SecondTimeDerivative("timestep", ppp),
            SecondTimeDerivative("Norm_Pressure", ppp),
            TimeIntegral("t", ppp),
            TimeIntegral("timestep", ppp),
            TimeIntegral("Norm_Pressure", ppp),
            ])

    # Attach a callback to postprocessor so we can inspect direct compute requests
    def ppcallback(field, data, t, timestep):
        print "DEBUG: In ppcallback", field.name, data, t, timestep
        ppcallback.calls[field.name].append((t, timestep))
    ppcallback.calls = defaultdict(list)
    pp._callback = ppcallback

    # Setup some mock problem state
    problem = MockProblem(dict(T=0.3, dt=0.1))
    dt = problem.params.dt
    T0 = problem.params.T0
    T = problem.params.T

    # Toy state evolution
    uexpr = Expression(("1.0 + 2*x[0] + 3*x[1]", "1.0 + 4*x[0] + 5*x[1]"))
    pexpr = Expression("t + 2.0*t*t + 3*t*t*t", t=0.0)
    #dpdtexpr = Expression("1.0 + 4.0*t + 9.0*t*t", t=0.0) # exact derivative
    #dpdtexpr = Expression("1.0 + 2.0*(2.0*t + dt) + 3.0*(3.0*t*t + 3.0*t*dt + dt*dt)", t=0.0, dt=0.0) # approximate derivative
    dpdtexpr = Expression("1.0 + 2.0*(2.0*(t-dt) + dt) + 3.0*(3.0*(t-dt)*(t-dt) + 3.0*(t-dt)*dt + dt*dt)", t=0.0, dt=0.0) # approximate derivative
    dpdtexpr.dt = dt

    # Setup some mock scheme state
    spaces = NSSpacePoolSplit(problem.mesh, 1, 1)
    u = Function(spaces.V)
    p = Function(spaces.Q)

    # Set "initial condition"
    u.interpolate(uexpr)

    # Update postprocessor for a number of timesteps, this is where the main code under test is
    for (t, timestep) in [(T0+i*dt, i) for i in range(int(0.5+T/dt)+1)]:
        # Fake a varying pressure
        pexpr.t = t
        p.interpolate(pexpr)
        assert abs(assemble(p**2*dx(domain=problem.mesh)) - (t+2*t**2+3*t**3)**2)  <  1e-8
        #assert assemble(dpdtexpr**2*dx, mesh=problem.mesh) < 1e-8

        # Run postprocessing step
        pp.update_all({"Velocity": u, "Pressure": p}, t, timestep, spaces, problem)

        # Test value of pressure derivative
        pr = pp.get("Pressure")
        #print "ZZZ", assemble(pr**2*dx, mesh=problem.mesh)
        if timestep > 0:
            dpdtexpr.t = t
            prm = pp.get("Pressure", -1)
            tdp = pp.get("TimeDerivative_Pressure")
            #print "IDS", pr, prm
            #print "YYY", assemble(prm**2*dx, mesh=problem.mesh)
            #print "YYY", assemble(dpdtexpr**2*dx, mesh=problem.mesh)
            #print "YYY", assemble(tdp**2*dx, mesh=problem.mesh)
            diff = assemble((dpdtexpr-tdp)**2*dx(domain=problem.mesh))
            assert abs(diff) < 1e-8

    # TODO: Check that we get the right amount of calls:
    #assert ppcallback.calls == defaultdict(list)

    # Get and check values from the final timestep
    assert abs( (pp.get("TimeDerivative_t")) - (1.0) ) < 1e-8
    assert abs( (pp.get("TimeIntegral_t")) - (0.5*T**2) ) < 1e-8
    assert abs( (pp.get("SecondTimeDerivative_t")) - (0.0) ) < 1e-8

    assert abs( (pp.get("TimeDerivative_timestep")) - (1.0/dt) ) < 1e-8
    # ts = (t-T0)/dt
   
    assert abs( (pp.get("SecondTimeDerivative_timestep")) - (0.0) ) < 1e-8

    # FIXME:
    #assert abs( (pp.get("TimeIntegral_timestep")) - (T/dt) ) < 1e-8
    #assert abs( (pp.get("TimeDerivative_Norm_Pressure")) - (0.0) ) < 1e-8
    #assert abs( (pp.get("TimeIntegral_Norm_Pressure")) - (0.0) ) < 1e-8
    #assert abs( (pp.get("SecondTimeDerivative_Norm_Pressure")) - (0.0) ) < 1e-8

def test_get_second_time_derivative():
    assert 1 == 1 # FIXME

def test_get_time_integral():
    assert 1 == 1 # FIXME

def test_get_norms():
    assert 1 == 1 # FIXME
