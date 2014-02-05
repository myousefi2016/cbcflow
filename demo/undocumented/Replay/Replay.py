from os import path
import sys

# Add FlowAroundCylinder problem as example problem
sys.path.insert(0, path.join(path.dirname(path.realpath(__file__)),'../../../demo/undocumented/FlowAroundCylinder'))

from cbcflow import *
from dolfin import *

from FlowAroundCylinder import FlowAroundCylinder       

def play():
    # First solve the problem
    problem = FlowAroundCylinder({"refinement_level": 3})
    scheme = IPCS_Stable()
    
    postprocessor = NSPostProcessor({"casedir": "Results"})
    
    postprocessor.add_fields([
        Velocity({"save": True, "stride_timestep": 2, "plot": True}),
        Pressure({"save": True, "stride_timestep": 3}),
    ])
    
    solver = NSSolver(problem, scheme, postprocessor)
    solver.solve()
    
def replay():
    # Create postprocessor pointing to the same casedir
    postprocessor = NSPostProcessor({"casedir": "Results"})
    
    # Add new fields to compute
    postprocessor.add_fields([
        StreamFunction({"save": True, "plot": True}),
        L2norm("Velocity", {"save": True, "plot": True}),
    ])
    
    # Replay
    replayer = NSReplay(postprocessor)
    replayer.replay()
    
if __name__ == '__main__':
    set_log_level(100)
    
    play()
    replay()
    