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
__author__ = "Oeyvind Evju <oyvinev@simula.no> and Martin Alnaes <martinal@simula.no>"
__date__ = "2013-04-26"
__copyright__ = "Copyright (C) 2013-2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from ..dol import plot, parameters, as_vector, MPI, Mesh

from time import time
import os
import glob
import re

from .paramdict import ParamDict
from .parameterized import Parameterized
from .utils import get_memory_usage, time_to_string, cbcflow_print, cbcflow_warning
from .restart import Restart

class NSSolver(Parameterized):
    "High level Navier-Stokes solver."
    def __init__(self, problem, scheme=None, postprocessor=None, params=None):
        Parameterized.__init__(self, params)

        self.problem = problem
        self.scheme = scheme
        self.postprocessor = postprocessor

    @classmethod
    def default_params(cls):
        # TODO: Insert generic nssolver params here
        params = ParamDict(
            debug=False,
            plot_solution=False,
            check_mem_frequency=0,
            restart = False, # Not yet supported
            restart_time = -1.0, # Not yet supported
            restart_timestep = -1, # Not yet supported
            enable_annotation=False,
            )
        return params

    def solve(self):
        self._reset()

        params = ParamDict(solver=self.params,
                           problem=self.problem.params,
                           scheme=self.scheme.params,
                           postprocessor=self.postprocessor.params)
        self.postprocessor.store_params(params)
        assert hasattr(self.problem, "mesh") and isinstance(self.problem.mesh, Mesh), "Unable to find problem.mesh!"
        self.postprocessor.store_mesh(self.problem.mesh)

        if self.params.restart:
            Restart(self.problem, self.postprocessor, self.params.restart_time, self.params.restart_timestep)

        # FIXME: Handle restart stuff, get from old ns script (see scratch/oldscripts/ns)
        # FIXME: Pick other details to reuse from old ns script, here or other places

        scheme_namespace = self.scheme.solve(self.problem, self.update)

        spaces = scheme_namespace["spaces"]
        self.postprocessor.finalize_all(spaces, self.problem)

        self._summarize()
        
        return scheme_namespace

    def _reset(self):
        self._initial_time = time()
        self._time = time()
        self._accumulated_time = 0

        if self.params.check_mem_frequency > 0:
            self._initial_memory = get_memory_usage()

    def _summarize(self):
        final_time = time()
        msg = "Total time spent in NSSolver: %s" % time_to_string(final_time - self._initial_time)
        cbcflow_print(msg)

        if self.params.check_mem_frequency > 0:
            final_memory = get_memory_usage()
            msg = "Memory usage before solve: %s\nMemory usage after solve: %s" % (
                self._initial_memory, final_memory)
            cbcflow_print(msg)

    def _update_timing(self, timestep, t, time_at_top):
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
        fr = self.params.check_mem_frequency
        if fr > 0 and timestep % fr == 0:
            # TODO: Report to file separately for each process
            cbcflow_print('Memory usage is: %s' % get_memory_usage())

    def _update_plot(self, u, p):
        if self.params.plot_solution:
            if MPI.num_processes() > 1:
                cbcflow_print("Unable to plot in parallel.")
                return

            if self.scheme._segregated:
                # TODO: If postprocessor makes a vector-valued u, we can use that here instead
                # TODO: Fix such that as_vector works for segregated schemes
                # TODO: Keep velocity plot in single window for segregated schemes (u not Variable, update not working)
                plot(u, title="Velocity")
            else:
                if not hasattr(self, 'uplot'):
                    self.uplot = plot(u, title="Velocity")
                else:
                    self.uplot.plot(u)
            if not hasattr(self, 'pplot'):
                self.pplot = plot(p, title="Pressure")
            else:
                self.pplot.plot(p)

    def update(self, u, p, t, timestep, spaces):
        # Do not run this if restarted from this timestep
        if self.params.restart and timestep==problem.params.start_timestep:
            return

        # Make a record of when update was called
        time_at_top = time()

        # Disable dolfin-adjoint annotation during the postprocessing
        if self.params.enable_annotation:
            parameters["adjoint"]["stop_annotating"] = True

        # Run postprocessor
        if self.postprocessor:
            self.postprocessor.update_all({"Velocity": u, "Pressure": p}, t, timestep, spaces, self.problem)

        # Plot solution
        self._update_plot(u, p)

        # Check memory usage
        self._update_memory(timestep)

        # Update timing data
        self._update_timing(timestep, t, time_at_top)

        # Enable dolfin-adjoint annotation before getting back to solve
        if self.params.enable_annotation:
            parameters["adjoint"]["stop_annotating"] = False
