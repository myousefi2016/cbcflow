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
from ..bases.MetaPPField import MetaPPField
from dolfin import *

class TimeIntegral(MetaPPField):
    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)

        t1 = pp.get("t")
        t0 = pp.get("t", -1)
        dt = t1 - t0

        if isinstance(u, Function):
            # Create placeholder for sum the first time
            if not hasattr(self, "_sum"):
                self._sum = Function(u.function_space())
            # Accumulate using backward euler integration
            self._sum.vector().axpy(dt, u.vector()) # FIXME: Validate this, not tested!
        else:
            # Create placeholder for sum the first time
            if not hasattr(self, "_sum"):
                self._sum = 0.0
            # Accumulate using backward euler integration
            self._sum += u

        # Increase count and return sum
        if not hasattr(self, "_count"):
            self._count = 0
        self._count += 1
        return self._sum
