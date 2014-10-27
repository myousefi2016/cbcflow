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

from cbcpost import SpacePool
from cbcpost import Field
from dolfin import Function

class PressureError(Field):
    @classmethod
    def default_params(cls):
        params = Field.default_params()
        params.replace(
            expr2function="project", # "assemble" | "project" | "interpolate"
            )
        return params

    def before_first_compute(self, get):
        p = get("Pressure")
        Q = p.function_space()
        spaces = SpacePool(Q.mesh())
        if self.params.expr2function == "assemble":
            degree = 0
        else:
            degree = Q.ufl_element().degree() + 1 # TODO: Is +1 sufficient?
        V = spaces.get_space(degree, 0)
        self._function = Function(V, name=self.name)

    def compute(self, get):
        p = get("Pressure")
        t = get("t")

        ua, pa = problem.analytical_solution(spaces, t)
        pe = pa - p

        return self.expr2function(pe, self._function)
