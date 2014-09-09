.. _Replay:

Replay a problem
=======================================

This demo is based on the :ref:`FlowAroundCylinder` demo, and demonstrates how one
can compute any :class:`Field` from a stored solution, given that the dependencies
of the field are either saved to disk or computable from the fields saved to disk.

This can be very useful when one needs additional information to the one specified
at the time of the solve.

The source code for this demo can be found in :download:`Replay.py`.

We start by importing from :ref:`FlowAroundCylinder`: ::

    import sys
    # Add FlowAroundCylinder problem as example problem
    sys.path.insert(0, path.join(path.dirname(path.realpath(__file__)),'../../../demo/documented/FlowAroundCylinder'))
    from FlowAroundCylinder import FlowAroundCylinder
    
We then import cbcflow ::

    from cbcflow import *
    
Note that we for this demo does not need to import from dolfin.

Play
____________________________________

We start by defining the *play*-routine, that is, a normal solve.
We use the :class:`.IPCS_Stable` scheme and save the results to *Results*. ::

    def play():
        # First solve the problem
        problem = FlowAroundCylinder({"refinement_level": 3})
        scheme = IPCS_Stable()


For this we simply store the velocity at every second timestep, and the
pressure at every third timestep. We plot the velocity, and we send keyword
*mode=color* so the plot will show the velocity magnitude: ::

        
        postprocessor = PostProcessor({"casedir": "Results"})
        
        postprocessor.add_fields([
            SolutionField("Velocity", {"save": True, "stride_timestep": 2, "plot": True, "plot_args": {"mode": "color"}}),
            SolutionField("Pressure", {"save": True, "stride_timestep": 3}),
        ])
        
We then solve the problem: ::

        solver = NSSolver(problem, scheme, postprocessor)
        solver.solve()


Replay
____________________________________
When the solve has been performed, we might want to compute some more derived fields.
We start by creating a new postprocessor instance. This instance must point to the same
case directory as the original solve: ::
    
    def replay():
        # Create postprocessor pointing to the same casedir
        postprocessor = PostProcessor({"casedir": "Results"})

We then define some new fields to store, and adds them to the postprocessor ::
    
    # Add new fields to compute
    postprocessor.add_fields([
        Stress({"save": True}),
        StreamFunction({"save": True, "plot": True}),
        Norm("Velocity", {"save": True, "plot": True}),
    ])
    
The :class:`.StreamFunction` and :class:`.Norm` of :class:`.Velocity` only depend on the
velocity, and we expect this to be computed at the same timesteps we computed the velocity
in the original solve.

The :class:`.Stress`, however, depends on both the :class:`.Velocity` and the :class:`Pressure`.
Since the velocity was saved at every second timestep, and the pressure was only saved every
third timestep, we expect this to be computed at timesteps 0, 6, 12, ... .

We then initiate a :class:`.NSReplay` instance with the postprocessor-instance as argument,
and call its *replay*-function to execute the replay routine: ::

    # Replay
    replayer = NSReplay(postprocessor)
    replayer.replay()

