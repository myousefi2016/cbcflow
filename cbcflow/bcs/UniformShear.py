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
from cbcflow.dol import *

from cbcflow.utils.bcs import compute_area

def compute_uniform_shear_value(u, ind, facet_domains, C=10000):
    mesh = facet_domains.mesh()
    A = compute_area(mesh, ind, facet_domains)
    dsi = ds[fd](ind)
    n = FacetNormal(mesh)
    d = len(n)
    u = as_vector(u)
    form = dot(u,n)*dsi
    Q = assemble(form)
    value = C*Q/A**1.5
    return value

# TODO: [martin] don't understand this design, why subclassing Constant?
class UniformShear(Constant):
    def __init__(self, u, ind, facet_domains, C=10000):
        Constant.__init__(self, compute_uniform_shear_value(u, ind, facet_domains, C))
