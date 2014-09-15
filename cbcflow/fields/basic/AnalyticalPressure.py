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
from cbcpost import Field
from dolfin import Function

class AnalyticalPressure(Field):
    @classmethod
    def default_params(cls):
        params = Field.default_params()
        params.replace(
            assemble=False,
            project=True,
            interpolate=False,
            )
        # TODO: Perhaps we should require that analytical_solution returns an Expression or Function and use interpolate instead?
        return params

    def before_first_compute(self, pp, spaces, problem):
        if self.params.assemble:
            degree = 0
        else:
            degree = spaces.p_degree + 1 # TODO: Is +1 sufficient?
        V = spaces.get_space(degree, 0)
        self._function = Function(V, name=self.name)

    def compute(self, pp, spaces, problem):
        t = pp.get("t")

        ua, pa = problem.analytical_solution(spaces, t)

        return self.expr2function(pa, self._function)
