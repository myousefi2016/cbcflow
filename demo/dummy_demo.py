#!/usr/bin/env python
# Hack to run without installing, useful while working
import sys
sys.path.insert(0,"../site-packages")

from headflow import *
from headflow.dol import *
set_log_level(100)
from numpy import linspace

class DummyProblem(NSProblem):
    def __init__(self):
        NSProblem.__init__(self)
        
        N = 8
        self.mesh = UnitCubeMesh(N, N, N)
        self.scale  = 2*(self.mesh.coordinates() - 0.5)
        self.mesh.coordinates()[:, :] = self.scale
        
        self.T = 1.0
        
        

class DummyScheme(NSScheme):
    def __init__(self):
        NSScheme.__init__(self)
        
        
        
    
    def solve(self, problem, update):
        V = VectorFunctionSpace(problem.mesh, "CG", 1)
        Q = FunctionSpace(problem.mesh, "CG", 1)
        
        u = Function(V)
        p = Function(Q)        
        
        # Velocity
        u_expr = Expression(('-((a*(pow(E,a*x[2])*cos(a*x[0] + d*x[1]) + pow(E,a*x[0])*sin(a*x[1] + d*x[2])))/pow(E,pow(d,2)*t*etabyrho))',
             '-((a*(pow(E,a*x[0])*cos(a*x[1] + d*x[2]) + pow(E,a*x[1])*sin(d*x[0] + a*x[2])))/pow(E,pow(d,2)*t*etabyrho))',
             '-((a*(pow(E,a*x[1])*cos(d*x[0] + a*x[2]) + pow(E,a*x[2])*sin(a*x[0] + d*x[1])))/pow(E,pow(d,2)*t*etabyrho))'),
            a=pi/4.0, d=2.0, E=e, rho=1.0, etabyrho=1.0, t=0.0)
        # Pressure
        p_expr = Expression('-(rho/2.0)*(pow(a,2)*(pow(E,2*a*x[0]) + pow(E,2*a*x[1]) + pow(E,2*a*x[2]) + 2*pow(E,a*(x[1] + x[2]))*cos(d*x[0] + a*x[2])*sin(a*x[0] + d*x[1]) + 2*pow(E,a*(x[0] + x[1]))*cos(a*x[1] + d*x[2])*sin(d*x[0] + a*x[2]) + 2*pow(E,a*(x[0] + x[2]))*cos(a*x[0] + d*x[1])*sin(a*x[1] + d*x[2])))/(pow(E,pow(d,2)*t*etabyrho))', a=pi/4.0, d=2.0, E=e, rho=1.0, etabyrho=1.0, t=0.0)

        t_range = linspace(0,1,201)
        for timestep, t in enumerate(t_range):
            u_expr.t = t
            p_expr.t = t
            
            u.assign(project(u_expr, V))
            p.assign(project(p_expr, Q))
            
            plot(u)
            plot(p)
            
            update(u,p,t,timestep)
            
class DummyPostProcessor(PostProcessorBase):
    def __init__(self, **kwargs):
        PostProcessorBase.__init__(self, **kwargs)
        
    def print_data(self):
        for inst in self.list_all:
            print inst.get_data()
            
    def print_all_params(self):
        for inst in self.list_all:
            print inst.params
            
            
class MockA(PPFieldBase):
   
    def __init__(self, **kwargs):
        PPFieldBase.__init__(self, **kwargs)
    
    def update(self, u, p, t, timestep, problem):
        value = timestep+10
        self.set_data(t, timestep, value)
        
        

class MockB(PPFieldBase):
    def __init__(self, **kwargs):
        assert("parent" in kwargs.keys())
        PPFieldBase.__init__(self, **kwargs)
       
    def update(self, u, p, t, timestep, problem):
        # Get parent data
        parent_datadict = self.parent.get_data()
        parent_data = parent_datadict["data"]
        
        value = parent_data**0.5
        
        value = {"sqrt": parent_data**0.5, "squared": parent_data**2, "cubed": parent_data**3}
        
        self.set_data(t, timestep, value)
        
    
class MockC(PPFieldBase):
    def __init__(self, **kwargs):
        PPFieldBase.__init__(self, **kwargs)

    def update(self, u, p, t, timestep, problem):
        # Get parent data
        value = u
        self.set_data(t, timestep, value)


        
if __name__ == '__main__':
    problem = DummyProblem()
    scheme = DummyScheme()
    postprocessor = DummyPostProcessor(casedir='dummy')
    
    timeparams = ParamDict(start_time=0.0, end_time=0.4, step_frequency=5)
    saveparams = ParamDict(save=True)

    
    a = MockA(timeparams=timeparams, saveparams=saveparams)
    b = MockB(parent=a, timeparams=timeparams, saveparams=saveparams)
    c = MockC(timeparams=timeparams, saveparams=saveparams)


    postprocessor.add_field(a)
    postprocessor.add_field(b)
    postprocessor.add_field(c)

    solver = NSSolver(problem, scheme, postprocessor)
    solver.solve()
    
    
    
            
