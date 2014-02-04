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
__author__ = "Oeyvind Evju <oyvinev@simula.no>"
__date__ = "2013-05-23"
__copyright__ = "Copyright (C) 2013-2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from .parameterized import Parameterized
from .paramdict import ParamDict
from cbcflow.dol import Constant

class TimeStepper(Parameterized): # TODO: Use this as base class for AdaptiveTimeStepper, see fixmes in that class
    def __init__(self, params):
        Parameterized.__init__(self, params)

        self.T = self.params.T
        self.T0 = self.params.T0
        self.dt = self.params.dt

        self.t = Constant(self.params.T0, name="TIME")
        self.timestep = self.params.timestep0

    @classmethod
    def default_params(cls):
        params = ParamDict(
            # Time parameters
            T = 1,
            dt = 1e-4,

            # Parameters for restart
            T0 = 0,
            timestep0 = 0,
            )
        return params

    def time(self):
        return self.t

    def current_timestep(self):
        return self.timestep

    def __iter__(self):
        return self

class ConstantTimestepping(TimeStepper): # TODO: Rename to ConstantTimeStepper
    def __init__(self, params):
        TimeStepper.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = TimeStepper.default_params()
        return params

    def next(self):
        eps = 1e-12
        if float(self.t) <= self.T - eps:
            self.t.assign(float(self.t) + self.dt)
            self.timestep += 1
            return self.timestep
        else:
            raise StopIteration

    def __len__(self):
        return int((self.T-self.T0)/self.dt)
