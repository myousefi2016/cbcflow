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
from cbcflow.fields.bases.MetaPPField import MetaPPField
from dolfin import Function

class TimeIntegral(MetaPPField):
    def compute(self, pp, spaces, problem):
        u1 = pp.get(self.valuename)
        u0 = pp.get(self.valuename, -1)

        t1 = pp.get("t")
        t0 = pp.get("t", -1)
        dt = t1 - t0

        if isinstance(u0, Function):
            # Create placeholder for sum the first time
            if not hasattr(self, "_sum"):
                self._sum = Function(u0.function_space())
            # Accumulate using trapezoidal integration
            self._sum.vector().axpy(dt/2.0, u0.vector()) # FIXME: Validate this, not tested!
            self._sum.vector().axpy(dt/2.0, u1.vector()) 
        else:
            # Create placeholder for sum the first time
            if not hasattr(self, "_sum"):
                self._sum = 0.0
            # Accumulate using trapezoidal integration
            self._sum += dt/2.0*u0
            self._sum += dt/2.0*u1

        return self._sum
