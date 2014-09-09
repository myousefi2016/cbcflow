from os import path
import sys
# Add Beltrami problem as example problem
sys.path.insert(0, path.join(path.dirname(path.realpath(__file__)), '../../../demo/undocumented/Beltrami'))
from Beltrami import Beltrami

from cbcflow import *



def play():
    problem = Beltrami(dict(dt=1e-2, T=1.0))
    scheme = IPCS()
    
    # Need to save velocity and pressure for restart to work
    fields = [
        SolutionField("Velocity", dict(save=True, stride_timestep=5)),
        SolutionField("Pressure", dict(save=True, stride_timestep=10)),
        Norm("Velocity", dict(save=True, stride_timestep=2))
    ]
    
    postprocessor = PostProcessor(dict(casedir='results'))
    
    postprocessor.add_fields(fields)
    
    solver = NSSolver(problem, scheme, postprocessor)
    solver.solve()


def restart():
    # Load params, to reuse
    import pickle
    params = pickle.load(open('results/params.pickle', 'r'))
    
    # Create new problem and scheme instances
    # Note: New scheme, and new end time
    problem = Beltrami(params.problem)
    problem.params.T = 2.0
    problem.params.dt = 5e-3
    scheme = IPCS_Stable(params.scheme)
    
    # Set up postprocessor with new fields
    fields = [
        SolutionField("Velocity", dict(save=True)),
        SolutionField("Pressure", dict(save=True)),
        WSS(problem, dict(save=True)),
    ]
    postprocessor = PostProcessor(dict(casedir='results'))
    postprocessor.add_fields(fields)
    
    # Set restart cbcflow-data
    solver = NSSolver(problem, scheme, postprocessor)
    solver.params["restart"] = True
    solver.params["restart_time"] = 0.5

    solver.solve()
    
if __name__ == '__main__':
    play()
    restart()