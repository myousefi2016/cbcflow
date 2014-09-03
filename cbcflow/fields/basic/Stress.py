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
from cbcflow.utils.common import sigma
from dolfin import Function, Constant

class Stress(PPField):
    def before_first_compute(self, pp, spaces, problem):
        if self.params.assemble:
            V = spaces.get_space(0, 2)
        else:
            V = spaces.DV
        self._function = Function(V, name=self.name)

    def compute(self, pp, spaces, problem):
        u = pp.get("Velocity")
        p = pp.get("Pressure")

        mu = Constant(problem.dynamic_viscosity())

        expr = sigma(u, p, mu)
        #u*epsilon(u) - p*Identity(u.cell().d) # TODO: is this with negative pressure?

        return self.expr2function(expr, self._function)
