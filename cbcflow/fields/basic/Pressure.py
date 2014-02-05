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
from cbcflow.fields.bases.PPField import PPField
from dolfin import Function, FunctionAssigner
from cbcflow.utils.core import NSSpacePoolMixed, NSSpacePoolSplit, NSSpacePoolSegregated

class Pressure(PPField):
    def convert(self, pp, spaces, problem):
        # Hack to get given p in whatever format it has
        p = super(Pressure, self).convert(pp, spaces, problem)

        if not isinstance(p, Function):
            if not hasattr(self, "_p"):
                self._p = Function(spaces.Q)
                assert isinstance(spaces, NSSpacePoolMixed)
                self._assigner = FunctionAssigner(spaces.Q, spaces.W.sub(1))

            # Hack: p is a Indexed(Coefficient()),
            # get the underlying mixed function
            w = p.operands()[0]
            self._assigner.assign(self._p, w.sub(1))

            p = self._p

        assert isinstance(p, Function)
        return p

SolverPressure = Pressure
