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
from cbcflow.post.fieldbases.Field import Field
from dolfin import as_vector, Function

class VelocityError(Field):
    @classmethod
    def default_params(cls):
        params = Field.default_params()
        params.replace(
            assemble=False,
            project=True,
            interpolate=False,
            )
        return params

    def before_first_compute(self, get):
        if self.params.assemble:
            degree = 0
        else:
            degree = spaces.u_degree + 1 # TODO: Is +1 sufficient?
        V = spaces.get_space(degree, 1)
        self._function = Function(V, name=self.name)

    def compute(self, get):
        u = get("Velocity")
        t = get("t")

        ua, pa = problem.analytical_solution(spaces, t)
        ua = as_vector(ua)
        ue = ua - u

        return self.expr2function(ue, self._function)