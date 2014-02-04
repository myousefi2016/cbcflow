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
from dolfin import as_vector, Function

class AnalyticalVelocity(PPField):
    @classmethod
    def default_params(cls):
        params = PPField.default_params()
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
            degree = spaces.u_degree + 1 # TODO: Is +1 sufficient?
        V = spaces.get_space(degree, 1)
        self._function = Function(V, name=self.name)

    def compute(self, pp, spaces, problem):
        u = pp.get("Velocity")
        t = pp.get("t")

        ua, pa = problem.analytical_solution(spaces, t)
        ua = as_vector(ua)

        return self.expr2function(ua, self._function)
