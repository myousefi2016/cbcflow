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

from cbcpost import SpacePool, get_grad_space, Field
from cbcflow.utils.common import sigma
from cbcflow.core.nsproblem import NSProblem
from dolfin import Function, Constant

class Stress(Field):
    def __init__(self, problem, params=None, name="default", label=None):
        Field.__init__(self, params, name, label)
        assert isinstance(problem, NSProblem)
        self.problem = problem

    def before_first_compute(self, get):
        u = get("Velocity")
        V = get_grad_space(u)

        #spaces = SpacePool(u.function_space().mesh())
        #V = spaces.get_grad_space(u.function_space())
        #if self.params.expr2function == "assemble":
        #    V = spaces.get_space(0, 2)
        #else:
        #    V = spaces.DV

        self._function = Function(V, name=self.name)

    def compute(self, get):
        u = get("Velocity")
        p = get("Pressure")

        # FIXME: Implement get("mu") instead
        mu = self.problem.params.mu
        if isinstance(mu, (float, int)):
            mu = Constant(mu)

        expr = sigma(u, p, mu)
        #u*epsilon(u) - p*Identity(u.cell().d) # TODO: is this with negative pressure?

        return self.expr2function(expr, self._function)
