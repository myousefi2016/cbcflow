#!/usr/bin/env py.test
"""
Tests of postprocessing framework in cbcflow.
"""

from collections import defaultdict

from cbcflow import (ParamDict, NSProblem, NSPostProcessor, NSScheme,
    PPField, Velocity, Pressure, VelocityGradient, Strain, Stress, WSS,
    TimeDerivative, SecondTimeDerivative, TimeIntegral, L2norm)
from cbcflow.fields import *

from cbcflow.utils.core import NSSpacePoolSplit
from cbcflow.utils.schemes import compute_regular_timesteps

import pytest

import dolfin
dolfin.set_log_level(40)
#from dolfin import (UnitSquareMesh, Function, Expression, norm, errornorm, assemble, dx,
#                    interpolate, plot)
from dolfin import *
from math import sqrt
from numpy.random import random
from numpy import linspace

# Avoid negative norms caused by instable tensor representation:
dolfin.parameters["form_compiler"]["representation"] = "quadrature"
dolfin.parameters["form_compiler"]["quadrature_degree"] = 1

class MockProblem(NSProblem):
    def __init__(self, D, params=None):
        NSProblem.__init__(self, params)
        self.D = D
        if D == 2:
            mesh = UnitSquareMesh(6,6)
        elif D == 3:
            mesh = UnitCubeMesh(3,3,3)
        self.initialize_geometry(mesh)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            T = 2.0,
            dt = 0.1,
            mu = 0.1,
            rho = 0.9,
            )
        return params

@pytest.fixture(scope="module", params=[2,3])
def problem(request):
    problem = MockProblem(request.param, dict(T=2.0))
    return problem

@pytest.fixture(scope="function", autouse=True)
def set_problem_params(request, problem):
    problem.params.dt = request.getfuncargvalue('dt')

@pytest.fixture(scope="module")
def spaces(problem):
    return NSSpacePoolSplit(problem.mesh, 1, 1)

@pytest.fixture(scope="function")
def pp():
    return NSPostProcessor()

class MockFunctionField(PPField):
    def before_first_compute(self, pp, spaces, problem):
        Q = spaces.Q
        t = pp.get('t')
        self.f = Function(Q)
        self.expr = Expression("1+x[0]*x[1]*t", t=t)
        
    def compute(self, pp, spaces, problem):
        t = pp.get('t')
        self.expr.t = t
        self.f.interpolate(self.expr)
        return self.f

class MockVectorFunctionField(PPField):
    def before_first_compute(self, pp, spaces, problem):
        V = spaces.V
        t = pp.get('t')
        self.f = Function(V)
        
        D = V.mesh().geometry().dim()
        if D == 2:
            self.expr = Expression(("1+x[0]*t", "3+x[1]*t"), t=t)
        elif D == 3:
            self.expr = Expression(("1+x[0]*t", "3+x[1]*t", "10+x[2]*t"), t=t)
        
    def compute(self, pp, spaces, problem):
        t = pp.get('t')
        self.expr.t = t
        self.f.interpolate(self.expr)
        return self.f
     
class MockTupleField(PPField):
    def compute(self, pp, spaces, problem):
        t = pp.get('t')
        return (t, 3*t, 1+5*t)

def test_TimeDerivative(problem, spaces, pp, start_time, end_time, dt):   
    # Setup some mock scheme state
    dt, timesteps, start_timestep = compute_regular_timesteps(problem)
    
    params = dict(finalize=True, start_time=start_time, end_time=end_time)
    
    pp.add_fields([
        MockFunctionField(),
        MockVectorFunctionField(),
        MockTupleField(),
    ])
    
    pp.add_fields([
            TimeDerivative("t", params),
            TimeDerivative("timestep", params),
            TimeDerivative("MockFunctionField", params),
            TimeDerivative("MockVectorFunctionField", params),
            TimeDerivative("MockTupleField", params),
            ])

    # Update postprocessor for a number of timesteps, this is where the main code under test is
    for timestep, t in enumerate(timesteps, start_timestep):
        # Run postprocessing step
        pp.update_all({}, t, timestep, spaces, problem)

    pp.finalize_all(spaces, problem)

    # Get and check values from the final timestep
    assert abs( (pp.get("TimeDerivative_t", compute=False)) - (1.0) ) < 1e-8
    assert abs( (pp.get("TimeDerivative_timestep", compute=False)) - (1.0/dt) ) < 1e-8
    assert errornorm(pp.get("TimeDerivative_MockFunctionField"), interpolate(Expression("x[0]*x[1]"), spaces.Q)) < 1e-8
    D = problem.D
    if D == 2:
        assert errornorm(pp.get("TimeDerivative_MockVectorFunctionField"), interpolate(Expression(("x[0]", "x[1]")), spaces.V)) < 1e-8
    elif D == 3:
        assert errornorm(pp.get("TimeDerivative_MockVectorFunctionField"), interpolate(Expression(("x[0]", "x[1]", "x[2]")), spaces.V)) < 1e-8
    
    assert max( [ abs(x1-x0) for x1,x0 in zip(pp.get("TimeDerivative_MockTupleField"), (1,3,5)) ] ) < 1e-8

def test_TimeIntegral(problem, spaces, pp, start_time, end_time, dt):
    # Setup some mock scheme state
    dt, timesteps, start_timestep = compute_regular_timesteps(problem)

    params = dict(finalize=True, start_time=start_time, end_time=end_time)
    
    pp.add_fields([
        MockFunctionField(),
        MockVectorFunctionField(),
        MockTupleField(),
    ])
    
    pp.add_fields([
            TimeIntegral("t", params),
            TimeIntegral("MockFunctionField", params),
            TimeIntegral("MockVectorFunctionField", params),
            TimeIntegral("MockTupleField", params),
            ])

    # Update postprocessor for a number of timesteps, this is where the main code under test is
    for timestep, t in enumerate(timesteps, start_timestep):
        # Run postprocessing step
        pp.update_all({}, t, timestep, spaces, problem)

    pp.finalize_all(spaces, problem)

    assert abs( pp.get("TimeIntegral_t") - (0.5*(end_time**2-start_time**2)) ) < 1e-8   
    assert errornorm(
        pp.get("TimeIntegral_MockFunctionField"),
        interpolate(Expression("t1-t0+0.5*x[0]*x[1]*(t1*t1-t0*t0)", t1=end_time, t0=start_time), spaces.Q)
        ) < 1e-8
    
    D = problem.D
    if D == 2:
        assert errornorm(
           pp.get("TimeIntegral_MockVectorFunctionField"),
           interpolate(Expression(("t1-t0+0.5*x[0]*(t1*t1-t0*t0)", "3*(t1-t0)+0.5*x[1]*(t1*t1-t0*t0)"), t1=end_time, t0=start_time), spaces.V)
        ) < 1e-8
    elif D == 3:
        assert errornorm(
           pp.get("TimeIntegral_MockVectorFunctionField"),
           interpolate(Expression(("t1-t0+0.5*x[0]*(t1*t1-t0*t0)", "3*(t1-t0)+0.5*x[1]*(t1*t1-t0*t0)", "10*(t1-t0)+0.5*x[2]*(t1*t1-t0*t0)"), t1=end_time, t0=start_time), spaces.V)
        ) < 1e-8
    
    I = 0.5*(end_time**2-start_time**2)
    assert max( [ abs(x1-x0) for x1,x0 in zip(pp.get("TimeIntegral_MockTupleField"), (I, 3*I, end_time-start_time+5*I) ) ] ) < 1e-8

def test_TimeAverage(problem, spaces, pp, start_time, end_time, dt):
    # Setup some mock scheme state
    dt, timesteps, start_timestep = compute_regular_timesteps(problem)
    
    params = dict(start_time=start_time, end_time=end_time, finalize=True)
    pp.add_fields([
        MockFunctionField(),
        MockVectorFunctionField(),
        MockTupleField(),
    ])
    
    pp.add_fields([
            TimeAverage("t", params),
            TimeAverage("MockFunctionField", params),
            TimeAverage("MockVectorFunctionField", params),
            TimeAverage("MockTupleField", params),
            ])

    # Update postprocessor for a number of timesteps, this is where the main code under test is
    for timestep, t in enumerate(timesteps, start_timestep):
        # Run postprocessing step
        pp.update_all({}, t, timestep, spaces, problem)

    pp.finalize_all(spaces, problem)
    
    assert abs( pp.get("TimeAverage_t") - (0.5*(end_time**2-start_time**2))/(end_time-start_time) ) < 1e-8
    assert errornorm(
        pp.get("TimeAverage_MockFunctionField"),
        interpolate(Expression("1+0.5*x[0]*x[1]*(t1+t0)", t1=end_time, t0=start_time), spaces.Q)
        ) < 1e-8
    
    D = problem.D
    if D == 2:
        assert errornorm(
           pp.get("TimeAverage_MockVectorFunctionField"),
           interpolate(Expression(("1+0.5*x[0]*(t1+t0)", "3+0.5*x[1]*(t1+t0)"), t1=end_time, t0=start_time), spaces.V)
        ) < 1e-8
    elif D == 3:
        assert errornorm(
           pp.get("TimeAverage_MockVectorFunctionField"),
           interpolate(Expression(("1+0.5*x[0]*(t1+t0)", "3+0.5*x[1]*(t1+t0)", "10+0.5*x[2]*(t1+t0)"), t1=end_time, t0=start_time), spaces.V)
        ) < 1e-8
    
    I = (0.5*end_time*end_time-0.5*start_time*start_time)/(end_time-start_time)    
    assert max( [ abs(x1-x0) for x1,x0 in zip(pp.get("TimeAverage_MockTupleField"), (I, 3*I, 1+5*I) ) ] ) < 1e-8
    
def test_TimeIntegral_of_TimeDerivative(problem, spaces, pp, start_time, end_time, dt):
    # Setup some mock scheme state
    dt, timesteps, start_timestep = compute_regular_timesteps(problem)
    
    pp.add_fields([
        MockFunctionField(),
        MockVectorFunctionField(),
        MockTupleField(),
    ])

    pp.add_fields([
        TimeDerivative("t"),
        TimeDerivative("MockFunctionField"),
        TimeDerivative("MockVectorFunctionField"),
        TimeDerivative("MockTupleField"),
        ])

    params = dict(start_time=start_time, end_time=end_time, finalize=False)
    pp.add_fields([
        TimeIntegral("TimeDerivative_t", params),
        TimeIntegral("TimeDerivative_MockFunctionField", params),
        TimeIntegral("TimeDerivative_MockVectorFunctionField", params),
        TimeIntegral("TimeDerivative_MockTupleField", params),
    ])
    
    # Because of backward finite differencing in TimeDerivative, a numerical error will be introduced if start_time<dt
    err_factor = max(0.0, dt-start_time)

    err_t = err_factor*0.5*(1-start_time/dt)
    err_MockFunctionField = err_factor*norm(interpolate(Expression("0.5*x[0]*x[1]*(1-t0/dt)", dt=dt, t0=start_time), spaces.Q))
    
    D = problem.mesh.geometry().dim()
    if D == 2:
        err_MockVectorFunctionField = err_factor*norm(interpolate(Expression(("0.5*x[0]*(1-t0/dt)", "0.5*x[1]*(1-t0/dt)"), dt=dt, t0=start_time), spaces.V))
    elif D == 3:
        err_MockVectorFunctionField = err_factor*norm(interpolate(Expression(("0.5*x[0]*(1-t0/dt)", "0.5*x[1]*(1-t0/dt)", "0.5*x[2]*(1-t0/dt)"), dt=dt, t0=start_time), spaces.V))
    else:
        raise Exception("D must be 2 or 3")
    
    err_MockTupleField = tuple([err_factor*(1-start_time/dt)*x/2.0 for x in [1.0, 3.0, 5.0]])
    
    if start_time > dt:
        assert err_factor == 0.0
    
     # Update postprocessor for a number of timesteps, this is where the main code under test is
    for timestep, t in enumerate(timesteps, start_timestep):
        # Run postprocessing step
        pp.update_all({}, t, timestep, spaces, problem)
        if start_time < t < end_time: 
            assert abs( pp.get("TimeIntegral_TimeDerivative_t") - ((t-start_time)-err_t) ) < 1e-8
    
    pp.finalize_all(spaces, problem)

    assert err_t-1e-8 < abs( pp.get("TimeIntegral_TimeDerivative_t") - (end_time-start_time)) < err_t+1e-8    
    assert err_MockFunctionField-1e-8 < \
            errornorm(pp.get("TimeIntegral_TimeDerivative_MockFunctionField"),
                    interpolate(Expression("x[0]*x[1]*(t1-t0)", t0=start_time, t1=end_time), spaces.Q)
                   ) < \
          err_MockFunctionField+1e-8

    if D == 2:
        assert err_MockVectorFunctionField-1e-8 < \
            errornorm(pp.get("TimeIntegral_TimeDerivative_MockVectorFunctionField"),
                    interpolate(Expression(("x[0]*(t1-t0)", "x[1]*(t1-t0)"), t0=start_time, t1=end_time), spaces.V)
                   ) < \
          err_MockVectorFunctionField+1e-8
    elif D == 3:
        assert err_MockVectorFunctionField-1e-8 < \
            errornorm(pp.get("TimeIntegral_TimeDerivative_MockVectorFunctionField"),
                    interpolate(Expression(("x[0]*(t1-t0)", "x[1]*(t1-t0)", "x[2]*(t1-t0)"), t0=start_time, t1=end_time), spaces.V)
                   ) < \
          err_MockVectorFunctionField+1e-8
    else:
        raise Exception("D must be 2 or 3")
    
    I = end_time-start_time
    assert abs( abs(pp.get("TimeIntegral_TimeDerivative_MockTupleField")[0] - I) - err_MockTupleField[0])  < 1e-8
    assert abs( abs(pp.get("TimeIntegral_TimeDerivative_MockTupleField")[1] - 3*I) - err_MockTupleField[1])  < 1e-8
    assert abs( abs(pp.get("TimeIntegral_TimeDerivative_MockTupleField")[2] - 5*I) - err_MockTupleField[2])  < 1e-8

def test_Maximum(problem, spaces, pp, start_time, end_time, dt):
        # Setup some mock scheme state
    dt, timesteps, start_timestep = compute_regular_timesteps(problem)
    
    pp.add_fields([
        MockFunctionField(),
        MockVectorFunctionField(),
        MockTupleField(),
    ])

    pp.add_fields([
        Maximum("t"),
        Maximum("MockFunctionField"),
        Maximum("MockVectorFunctionField"),
        Maximum("MockTupleField"),
        ])
    
    xmax = MPI.max(max(problem.mesh.coordinates()[:,0]))
    ymax = MPI.max(max(problem.mesh.coordinates()[:,1]))
    if problem.D > 2:
        zmax = MPI.max(max(problem.mesh.coordinates()[:,2]))
    
     # Update postprocessor for a number of timesteps, this is where the main code under test is
    for timestep, t in enumerate(timesteps, start_timestep):
        # Run postprocessing step
        pp.update_all({}, t, timestep, spaces, problem)
        if start_time < t < end_time:
            assert abs(pp.get("Maximum_t") - t) < 1e-8
            assert abs(pp.get("Maximum_MockFunctionField") - (1+xmax*ymax*t)) < 1e-8
            if problem.D == 2:
                assert abs(pp.get("Maximum_MockVectorFunctionField") - (3+ymax*t)) < 1e-8
            elif problem.D == 3:
                assert abs(pp.get("Maximum_MockVectorFunctionField") - (10+zmax*t)) < 1e-8
            assert abs(pp.get("Maximum_MockTupleField") - (1+5*t)) < 1e-8
    
    pp.finalize_all(spaces, problem)
    
    assert abs(pp.get("Maximum_t") - t) < 1e-8
    assert abs(pp.get("Maximum_MockFunctionField") - (1+xmax*ymax*t)) < 1e-8
    if problem.D == 2:
        assert abs(pp.get("Maximum_MockVectorFunctionField") - (3+ymax*t)) < 1e-8
    elif problem.D == 3:
        assert abs(pp.get("Maximum_MockVectorFunctionField") - (10+zmax*t)) < 1e-8
    assert abs(pp.get("Maximum_MockTupleField") - (1+5*t)) < 1e-8
    
    
    



