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

from dolfin import Function, grad, det

class Delta(PPField):
    @classmethod
    def default_params(cls):
        params = PPField.default_params()
        params.replace(
            assemble=False, # Change to this to use assemble into DG0 by default
            project=True,
            interpolate=False,
            )
        return params

    def before_first_compute(self, pp, spaces, problem):
        if self.params.assemble:
            V = spaces.get_space(0, 0)
        else:
            # Accurate degree is 6*(spaces.u_degree-1)
            degree = 1
            V = spaces.get_space(degree, 0)
        self._function = Function(V, name=self.name)

    def compute(self, pp, spaces, problem):
        u = pp.get("Velocity")
        Q = pp.get("Q")

        #S = (grad(u) + grad(u).T)/2
        #Omega = (grad(u) - grad(u).T)/2
        expr = (Q**3 / 3 + det(grad(u))**2 / 2 )

        return self.expr2function(expr, self._function)
