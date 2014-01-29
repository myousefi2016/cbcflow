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

from ..bases.MetaPPField2 import MetaPPField2
from dolfin import *

class DiffL2norm(MetaPPField2):
    "Compute the L2 norm of the difference between uh and u relative to u."
    def compute(self, pp, spaces, problem):
        uh = pp.get(self.valuename1)
        u = pp.get(self.valuename2)
        e = uh - u

        dx = problem.dx

        uh2 = assemble(uh**2*dx(), mesh=problem.mesh)
        u2 = assemble(u**2*dx(), mesh=problem.mesh)
        e2 = assemble(e**2*dx(), mesh=problem.mesh)

        # Relative norm
        eps = 1e-14
        if abs(u2) > eps:
            return sqrt(e2 / u2)
        else:
            return sqrt(e2)
