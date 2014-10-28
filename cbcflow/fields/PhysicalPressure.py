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
from dolfin import Function, FunctionAssigner
from cbcflow import NSProblem

class PhysicalPressure(Field):
    "The physical pressure is the solver pressure scaled by density."
    def __init__(self, problem, params=None, name="default", label=None):
        Field.__init__(self, params, name, label)
        assert isinstance(problem, NSProblem)
        self.problem = problem

    def before_first_compute(self, get):
        p = get("Pressure")

        if isinstance(p, Function):
            Q = p.function_space()
        else:
            # p is part of a mixed function w=(u,p)
            w = p.operands()[0]
            W = w.function_space()
            W1 = W.sub(1)
            Q = W1.collapse() # TODO: This creates a new function space
            self._assigner = FunctionAssigner(Q, W1)

        self._p = Function(Q)

    def compute(self, get):
        # Get solver pressure
        p = get("Pressure")

        if isinstance(p, Function):
            self._p.assign(p)
        else:
            # Hack: get the underlying subfunction from mixed function w
            #       (w = Function(V*Q); p = Indexed(w)):
            self._assigner.assign(self._p, p.operands()[0].sub(1))

        # Scale by density
        rho = self.problem.params.rho
        pv = self._p.vector()
        pv *= rho
        return self._p
