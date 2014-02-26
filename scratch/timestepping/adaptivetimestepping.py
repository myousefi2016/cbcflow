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

from cbcflow.core.parameterized import Parameterized
from cbcflow.core.paramdict import ParamDict
from cbcflow.dol import Constant

# TODO: Make base class, there's a bit of shared code here
class AdaptiveTimestepping(Parameterized): # TODO: Better name is AdaptiveTimeStepper
    def __init__(self, params=None):
        Parameterized.__init__(self, params)

        self.T = self.params.T
        self.T0 = self.params.T0
        self.dt = self.params.dt

        self.t = Constant(self.params.T0, name="TIME")
        self.timestep = self.params.timestep0

        # Set inital values such that checks allow initial increase
        self.last_modified_timestep = self.timestep - 2*self.params.min_timesteps_between_increases
        self.previous_check_value = (self.params.max_value + self.params.min_value) // 2

    @classmethod
    def default_params(cls):
        #params = TimeStepper.default_params()
        #params.update(
        params = ParamDict(
            # Time parameters
            T = 1,
            dt = 1e-4,

            # Parameters for restart
            T0 = 0,
            timestep0 = 0,

            # Timestep adjustment factors
            decrease_factor = 2.0,
            increase_factor = 1.5,

            # Timestep check parameters
            min_value = 15,
            max_value = 40,
            min_timesteps_between_increases = 5,
            increase_on_positive_gradient = False,
            )
        return params

    def time(self):
        return self.t

    def current_timestep(self):
        return self.timestep

    def __iter__(self):
        return self

    def next(self):
        eps = 1e-12
        if float(self.t) <= self.T - eps:
            self.t.assign(float(self.t) + self.dt)
            self.timestep += 1
            return self.timestep
        else:
            raise StopIteration

    def adjusted_timestep(self, check_value):
        increased = False
        decreased = False

        old_dt = self.dt

        if check_value > self.params.max_value:
            # Reset to previous timestep
            self.timestep -= 1
            self.t.assign(float(self.t) - old_dt)

            # Find new smaller dt
            decreased = True
            self.dt = old_dt / self.params.decrease_factor

        elif check_value < self.params.min_value:
            if self.timestep - self.last_modified_timestep > self.params.min_timesteps_between_increases:
                # Never increase more often than specified (NB! also waits this long after a decrease!)
                increased = False
            elif check_value <= self.previous_check_value:
                # Always allow increase if check value is decreasing
                # (i.e. moving further outside acceptance range)
                # This always triggers on the first timestep to move outside the acceptance range,
                # but may not be true for subsequent timesteps which may stay outside the range and
                # move towards acceptable values.
                increased = True
            else:
                # Allow increase independent of gradient sign?
                # If this is false we don't increase if things are going
                # the right way, to avoid getting too large dt quickly.
                increased = self.params.increase_on_positive_gradient

            if increased:
                # Find a new larger dt
                self.dt = old_dt * self.params.increase_factor

        # Bookkeeping
        is_modified = (decreased or increased)
        if is_modified:
            self.last_modified_timestep = self.timestep
        self.previous_check_value = check_value

        return (is_modified, self.dt)
