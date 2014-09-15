.. _Restart:

Restart a problem
=====================================
This demo demonstrates the functionality of restarting a problem. This could be useful
for example if one wants to run the simulation for longer than anticipated on initial
solve, if one wants parts of the solution to have a different temporal resolution or if
one wants to completely change the boundary conditions of the problem at a certain point.

It is based on the Beltrami problem, a 3D problem with a known analytical solution.

The source code for this demo can be found in :download:`Restart.py`.

We start by importing the problem, ::

    from os import path
    import sys
    # Add Beltrami problem as example problem
    sys.path.insert(0, path.join(path.dirname(path.realpath(__file__)), '../../../demo/undocumented/Beltrami'))
    from Beltrami import Beltrami

and cbcflow: ::

    from cbcflow import *


Play
______________________________________
We define the normal solve routine, by first initiating a problem and scheme: ::

    def play():
        problem = Beltrami(dict(dt=1e-2, T=1.0))
        scheme = IPCS()

Note that we here use the :class:`IPCS` scheme.


We then creates some fields to save for this solve: ::
        
        fields = [
            Velocity(dict(save=True, stride_timestep=5)),
            Pressure(dict(save=True, stride_timestep=10)),
            L2norm("Velocity", dict(save=True, stride_timestep=2))
        ]

Note that we *must* save velocity and pressure for restart to work.

We then add the fields to a :class:`.PostProcessor` instance, ::

        postprocessor = PostProcessor(dict(casedir='results'))
        postprocessor.add_fields(fields)

and solves the equation: ::
        
        solver = NSSolver(problem, scheme, postprocessor)
        solver.solve()



Restart
_______________________________________
When we restart the problem, we wish to reuse most of the parameters of the original
solve. These are save by the postprocessor to the case directory, and can be *unpickled*: ::
        
    def restart():
        # Load params, to reuse
        import pickle
        params = pickle.load(open('results/params.pickle', 'r'))
        
If we don't change any of the parameters, we would basically be solving the exact
same problem. Thus, we change the end time and time step of the problem: ::
        
        problem = Beltrami(params.problem)
        problem.params.T = 2.0
        problem.params.dt = 5e-3
        
We are also free to change the scheme, so we change the scheme to :class:`IPCS_Stable`: ::
        
        scheme = IPCS_Stable(params.scheme)
        
On the restart, we also set up a different set of fields ::

        # Set up postprocessor with new fields
        fields = [
            Velocity(dict(save=True)),
            Pressure(dict(save=True)),
            WSS(dict(save=True)),
        ]
        
and a new :class:`PostProcessor` instance: ::
        
        postprocessor = PostProcessor(dict(casedir='results'))
        postprocessor.add_fields(fields)
        
We then need to define our new solver, and set some restart-specific parameters: ::

    solver = NSSolver(problem, scheme, postprocessor)
    solver.params["restart"] = True
    solver.params["restart_time"] = 0.5

The solver will try to search for a solution in the postprocessors case directory at
time 0.5, and replace the the method :func:`.initial_conditions()` in the :class:`.NSProblem`
instance to reflect the solution at t=0.5.

Our call to solve will then restart this problem from the specified parameters,
and solve the problem: ::

    solver.solve()