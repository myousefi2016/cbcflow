from os import path
import sys

# Add Beltrami problem as example problem
sys.path.insert(0, path.join(path.dirname(path.realpath(__file__)), '../../../demo/undocumented/Beltrami'))

import pickle

from cbcflow import *
from dolfin import *

from Beltrami import Beltrami

def play():
    problem = Beltrami(dict(dt=1e-2, T=1.0))
    scheme = IPCS()
    
    # Need to save velocity and pressure for restart to work
    fields = [
        Velocity(dict(save=True, stride_timestep=5)),
        Pressure(dict(save=True, stride_timestep=10)),
        L2norm("Velocity", dict(save=True, stride_timestep=2))
    ]
    
    postprocessor = NSPostProcessor(dict(casedir='results'))
    
    postprocessor.add_fields(fields)
    
    solver = NSSolver(problem, scheme, postprocessor)
    solver.solve()


def restart():
    # Load params, to reuse
    params = pickle.load(open('results/params.pickle', 'r'))
    
    # Create new problem and scheme instances
    # Note: New scheme, and new end time
    problem = Beltrami(params.problem)
    problem.params.T = 2.0
    scheme = IPCS_Stable(params.scheme)
    
    # Set up postprocessor with new fields
    fields = [
        Velocity(dict(save=True)),
        Pressure(dict(save=True)),
        WSS(dict(save=True)),
    ]
    postprocessor = NSPostProcessor(dict(casedir='results'))
    postprocessor.add_fields(fields)
    
    # Set restart data
    solver_params = params.solver
    solver_params["restart"] = True
    solver_params["restart_time"] = 0.5
        
    solver = NSSolver(problem, scheme, postprocessor, params=solver_params)
    solver.solve()
    
if __name__ == '__main__':
    play()
    restart()