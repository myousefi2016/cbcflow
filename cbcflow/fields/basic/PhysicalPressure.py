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
from ..bases.PPField import PPField
from dolfin import Function
from ...core.spaces import NSSpacePoolMixed, NSSpacePoolSplit, NSSpacePoolSegregated

class PhysicalPressure(PPField):
    "The physical pressure is the solver pressure scaled by density."
    def compute(self, pp, spaces, problem):
        # Get solver pressure
        p = pp.get("Pressure")

        if not hasattr(self, "_p"):
            self._p = Function(spaces.Q)

        if isinstance(p, Function):
            self._p.assign(p)
        else:
            assert isinstance(spaces, NSSpacePoolMixed)
            if not hasattr(self, "_assigner"):
                self._assigner = FunctionAssigner(spaces.Q, spaces.W.sub(1))
            # Hack: p is a Indexed(Coefficient()),
            # get the underlying mixed function
            w = p.operands()[0]
            self._assigner.assign(self._p, w.sub(1))

        # Scale by density
        rho = problem.params.rho
        pv = self._p.vector()
        pv *= rho
        p = self._p

        assert isinstance(p, Function)
        return p
