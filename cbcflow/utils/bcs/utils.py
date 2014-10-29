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

from cbcflow.dol import (SubsetIterator, MPI, ds, assemble, Constant, sqrt,
                         FacetNormal, as_vector, mpi_comm_world, SpatialCoordinate)
import numpy as np

def x_to_r2(x, c, n):
    # TODO: Simplify this after testing
    d = len(c)
    rv = [x[i]-c[i] for i in xrange(d)]
    rvn = sum([rv[i]*n[i] for i in xrange(d)])
    rv = [rv[i] - rvn*n[i] for i in xrange(d)]
    r2 = sum(rv[i]**2 for i in xrange(d))
    return r2

def compute_radius(mesh, facet_domains, ind, center):
    d = len(center)
    it = SubsetIterator(facet_domains, ind)
    geom = mesh.geometry()
    #maxr2 = -1.0
    maxr2 = 0
    for i, facet in enumerate(it):
        ent = facet.entities(0)
        for v in ent:
            p = geom.point(v)
            r2 = sum((p[j] - center[j])**2 for j in xrange(d))
            maxr2 = max(maxr2, r2)
    r = MPI.max(mpi_comm_world(), sqrt(maxr2))
    return r

def compute_boundary_geometry_acrn(mesh, ind, facet_domains):
    # Some convenient variables
    assert facet_domains is not None
    dsi = ds(ind, domain=mesh, subdomain_data=facet_domains)

    d = mesh.geometry().dim()
    x = SpatialCoordinate(mesh)

    # Compute area of boundary tesselation by integrating 1.0 over all facets
    A = assemble(Constant(1.0)*dsi)
    assert A > 0.0, "Expecting positive area, probably mismatch between mesh and markers!"

    # Compute barycenter by integrating x components over all facets
    c = [assemble(x[i]*dsi) / A for i in xrange(d)]

    # Compute average normal (assuming boundary is actually flat)
    n = FacetNormal(mesh)
    ni = [assemble(n[i]*dsi) for i in xrange(d)]
    n_len = np.sqrt(sum([ni[i]**2 for i in xrange(d)])) # Should always be 1!?
    normal = [ni[i]/n_len for i in xrange(d)]

    # Compute radius by taking max radius of boundary points
    # (assuming boundary points are on exact geometry)
    r = compute_radius(mesh, facet_domains, ind, c)
    #r = np.sqrt(A / pi) # This old estimate is a few % lower because of boundary discretization errors

    return A, c, r, normal

def compute_area(mesh, ind, facet_domains):
    # Some convenient variables
    assert facet_domains is not None
    dsi = ds(ind, domain=mesh, subdomain_data=facet_domains)

    # Compute area of boundary tesselation by integrating 1.0 over all facets
    A = assemble(Constant(1.0)*dsi)
    assert A > 0.0, "Expecting positive area, probably mismatch between mesh and markers!"
    return A

def compute_transient_scale_value(bc, period, mesh, facet_domains, ind, scale_value):
    dsi = ds(ind, domain=mesh, subdomain_data=facet_domains)
    form = sqrt(as_vector(bc)**2) * dsi

    N = 100
    qt = [0]*N
    for i, t in enumerate(np.linspace(0, period, N)):
        for e in bc:
            e.set_t(t)
        qt[i] = assemble(form)
    for e in bc:
        e.set_t(0.0)

    q_avg = sum(qt) / len(qt)
    return scale_value / q_avg
