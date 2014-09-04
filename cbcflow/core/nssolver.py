# Copyright (C) 2010-2014 Simula Research Laboratory
#
# This file is part of CBCFLOW.
#
# CBCFLOW is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CBCFLOW is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CBCFLOW. If not, see <http://www.gnu.org/licenses/>.
from __future__ import division

from cbcflow.dol import parameters, Mesh, MPI

from time import time

from cbcflow.core.paramdict import ParamDict
from cbcflow.core.parameterized import Parameterized
from cbcflow.utils.common import get_memory_usage, time_to_string, cbcflow_print
from cbcflow.core.restart import Restart
from cbcflow.utils.common import Timer, PressureConverter, VelocityConverter

class NSSolver(Parameterized):
    """High level Navier-Stokes solver. This handles all logic between the cbcflow
    components.
    
    For full functionality, the user should instantiate this class with a NSProblem
    instance, NSScheme instance and NSPostProcessor instance.
    
    """
    def __init__(self, problem, scheme=None, postprocessor=None, params=None):
        Parameterized.__init__(self, params)

        self.problem = problem
        self.scheme = scheme
        self.postprocessor = postprocessor
        
        self.velocity_converter = VelocityConverter()
        self.pressure_converter = PressureConverter()

    @classmethod
    def default_params(cls):
        """Returns the default parameters for a problem.

        Explanation of parameters:
          - debug: bool, debug mode
          - check_memory_frequency: int, timestep frequency to check memory consumption
          - restart: bool, turn restart mode on or off
          - restart_time: float, time to search for restart data
          - restart_timestep: int, timestep to search for restart data
          
          If restart=True, maximum one of restart_time and restart_timestep can be set.
          
        """
        # TODO: Insert generic nssolver params here
        params = ParamDict(
            debug=False,
            check_memory_frequency=0,
            restart = False, 
            restart_time = -1.0,
            restart_timestep = -1,
            enable_annotation=False,
            timer_frequency=0,
            )
        return params

    def solve(self):
        """ Handles top level logic related to solve.
        
        Cleans casedir or loads restart data, stores parameters and mesh in
        casedir, calls scheme.solve, and lets postprocessor finalize all
        fields.
        
        Returns: namespace dict returned from scheme.solve
        """
        params = ParamDict(solver=self.params,
                           problem=self.problem.params,
                           scheme=self.scheme.params)
        
        
        self.timer = Timer(self.params.timer_frequency)
        self.timer._N = -1
        
        self._reset()
        
        if self.postprocessor != None:
            self.postprocessor._timer = self.timer
        
            if self.params.restart:
                Restart(self.problem, self.postprocessor, self.params.restart_time, self.params.restart_timestep)
                self.timer.completed("set up restart")
            else:
                # If no restart, remove any existing data coming from cbcflow
                self.postprocessor._clean_casedir()
                self.timer.completed("cleaned casedir")
            
            params = ParamDict(solver=self.params,
                               problem=self.problem.params,
                               scheme=self.scheme.params,
                               postprocessor=self.postprocessor.params)
            self.postprocessor.store_params(params)
            assert hasattr(self.problem, "mesh") and isinstance(self.problem.mesh, Mesh), "Unable to find problem.mesh!"
            self.postprocessor.store_mesh(self.problem.mesh)
            self.timer.completed("stored mesh")
            
        # FIXME: Pick other details to reuse from old ns script, here or other places

        scheme_namespace = self.scheme.solve(self.problem, self.update, self.timer)

        spaces = scheme_namespace["spaces"]
        
        if self.postprocessor != None:
            self.postprocessor.finalize_all()

        self._summarize()
        
        return scheme_namespace

    def _reset(self):
        "Reset timers and memory usage"
        self._initial_time = time()
        self._time = time()
        self._accumulated_time = 0

        if self.params.check_memory_frequency > 0:
            self._initial_memory = get_memory_usage()

    def _summarize(self):
        "Summarize time spent and memory (if requested)"
        final_time = time()
        msg = "Total time spent in NSSolver: %s" % time_to_string(final_time - self._initial_time)
        cbcflow_print(msg)

        if self.params.check_memory_frequency > 0:
            final_memory = get_memory_usage()
            msg = "Memory usage before solve: %s\nMemory usage after solve: %s" % (
                self._initial_memory, final_memory)
            cbcflow_print(msg)

    def _update_timing(self, timestep, t, time_at_top):
        "Update timing and print solver status to screen."
        # Time since last update equals the time for scheme solve
        solve_time = time_at_top - self._time

        # Store time for next update
        self._time = time()

        # Time since time at top of update equals postproc + some update overhead
        pp_time = self._time - time_at_top

        # Accumulate time spent for time left estimation
        if timestep > 1: # TODO: Call update for initial conditions, and skip timestep 0 instead for better estimates
            # (skip first step where jit compilation dominates)
            self._accumulated_time += solve_time
            self._accumulated_time += pp_time
        if timestep == 2:
            # (count second step twice to compensate for skipping first) (this may be overkill)
            self._accumulated_time += solve_time
            self._accumulated_time += pp_time

        # Report timing of this step
        if t:
            spent = time_to_string(self._accumulated_time)
            remaining = time_to_string(self._accumulated_time*(self.problem.params.T-t)/t)
        else:
            spent = '--'
            remaining = '--'
        msg = ("Timestep %5d finished (t=%.2e, %.1f%%) in %3.1fs (solve: %3.1fs). Time spent: %s Time remaining: %s" \
                                       % (timestep, t,
                                          100*t/self.problem.params.T,
                                          solve_time + pp_time,
                                          solve_time,
                                          spent,
                                          remaining,
                                          ))
        # TODO: Report to file, with additional info like memory usage, and make reporting configurable
        cbcflow_print(msg)

    def _update_memory(self, timestep):
        fr = self.params.check_memory_frequency
        if fr > 0 and timestep % fr == 0:
            # TODO: Report to file separately for each process
            cbcflow_print('Memory usage is: %s' % MPI.sum(get_memory_usage()))       

    def update(self, u, p, t, timestep, spaces):
        """Callback from scheme.solve after each timestep to handle update of
        postprocessor, timings, memory etc."""
        self.timer.completed("completed solve")
        
        # Do not run this if restarted from this timestep
        if self.params.restart and timestep==self.problem.params.start_timestep:
            return

        # Make a record of when update was called
        time_at_top = time()

        # Disable dolfin-adjoint annotation during the postprocessing
        if self.params.enable_annotation:
            parameters["adjoint"]["stop_annotating"] = True

        # Run postprocessor
        if self.postprocessor:
            #self.postprocessor.update_all({"Velocity": u, "Pressure": p}, t, timestep, spaces, self.problem)
            #self.postprocessor.update_all({"Velocity": u, "Pressure": p}, t, timestep)
            self.postprocessor.update_all({"Velocity": lambda: self.velocity_converter(u, spaces),
                                           "Pressure": lambda: self.pressure_converter(p, spaces)},
                                            t, timestep)
        self.timer.completed("postprocessor update")

        # Check memory usage
        self._update_memory(timestep)
        
        # Update timing data
        self.timer.increment()
        self._update_timing(timestep, t, time_at_top)

        # Enable dolfin-adjoint annotation before getting back to solve
        if self.params.enable_annotation:
            parameters["adjoint"]["stop_annotating"] = False
