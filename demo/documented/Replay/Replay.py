import sys
from os import path
# Add FlowAroundCylinder problem as example problem
sys.path.insert(0, path.join(path.dirname(path.realpath(__file__)),'../../../demo/documented/FlowAroundCylinder'))
from FlowAroundCylinder import FlowAroundCylinder       

from cbcflow import *

def play():
    # First solve the problem
    problem = FlowAroundCylinder({"refinement_level": 3})
    scheme = IPCS_Stable()
    
    postprocessor = PostProcessor({"casedir": "Results"})
    
    postprocessor.add_fields([
        Velocity({"save": True, "stride_timestep": 2, "plot": True, "plot_args": {"mode": "color"}}),
        Pressure({"save": True, "stride_timestep": 3}),
    ])
    
    solver = NSSolver(problem, scheme, postprocessor)
    solver.solve()
    
def replay():
    # Create postprocessor pointing to the same casedir
    postprocessor = PostProcessor({"casedir": "Results"})
    
    # Add new fields to compute
    postprocessor.add_fields([
        Stress({"save": True}),
        StreamFunction({"save": True, "plot": True}),
        Norm("Velocity", {"save": True, "plot": True}),
    ])
    
    # Replay
    replayer = NSReplay(postprocessor)
    replayer.replay()
    
if __name__ == '__main__':
    # Solve problem
    play()
    
    # Loop through saved solution and do some more calculations
    replay()
    