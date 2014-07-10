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
from dolfin import Function, FunctionAssigner, error
from cbcflow.utils.core import NSSpacePoolMixed, NSSpacePoolSegregated

class Velocity(PPField):
    def convert(self, pp, spaces, problem):
        # Hack to get given u in whatever format it has,
        # avoiding circular reference to this field      
        u = super(Velocity, self).convert(pp, spaces, problem)
        if u == None:
            return None
        d = spaces.d

        if not isinstance(u, Function):
            if not hasattr(self, "_u"):
                self._u = Function(spaces.V)

                if isinstance(spaces, NSSpacePoolMixed):
                    self._assigner = FunctionAssigner(spaces.V, spaces.W.sub(0))
                elif isinstance(spaces, NSSpacePoolSegregated):
                    self._assigner = FunctionAssigner(spaces.V, [spaces.U]*d)
                else:
                    error("It doesnt make sense to create a function assigner for a split space.")

            if isinstance(spaces, NSSpacePoolMixed):
                # Hack: u is a ListTensor([Indexed(Coefficient()),...]),
                # get the underlying mixed function
                w = u.operands()[0].operands()[0]
                assert w.shape() == (d+1,)
                us = w.sub(0)

            elif isinstance(spaces, NSSpacePoolSegregated):
                us = [u[i] for i in range(d)]

            self._assigner.assign(self._u, us)
            u = self._u

        assert isinstance(u, Function)
        return u
