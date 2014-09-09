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
from math import sqrt
from dolfin import assemble

class KineticEnergy(Field):
    def compute(self, get):
        u = get("Velocity")
        dx = problem.dx
        u_norms = [assemble(u[d]**2*dx()) for d in range(u.shape()[0])]
        energy = sqrt(sum(u_norms[d] for d in range(u.shape()[0])))
        return energy