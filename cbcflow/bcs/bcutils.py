from ..dol import *
import numpy as np

def x_to_r2(x, c, n):
    # TODO: Simplify this after testing
    d = len(c)
    rv = [x[i]-c[i] for i in xrange(d)]
    rvn = sum(rv[i]*n[i] for i in xrange(d))
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
    r = MPI.max(sqrt(maxr2))
    return r

def compute_boundary_geometry_acrn(mesh, ind, facet_domains):
    # Some convenient variables
    assert facet_domains is not None
    dsi = ds[facet_domains](ind)
    cell = mesh.ufl_cell()
    d = cell.d
    x = cell.x

    # Compute area of boundary tesselation by integrating 1.0 over all facets
    A = assemble(Constant(1.0)*dsi, mesh=mesh)
    assert A > 0.0, "Expecting positive area, probably mismatch between mesh and markers!"

    # Compute barycenter by integrating x components over all facets
    c = [assemble(x[i]*dsi, mesh=mesh) / A for i in xrange(d)]

    # Compute average normal (assuming boundary is actually flat)
    n = FacetNormal(mesh)
    ni = [assemble(n[i]*dsi, mesh=mesh) for i in xrange(d)]
    n_len = np.sqrt(sum(ni[i]**2 for i in xrange(d))) # Should always be 1!?
    normal = [ni[i]/n_len for i in xrange(d)]

    # Compute radius by taking max radius of boundary points
    # (assuming boundary points are on exact geometry)
    r = compute_radius(mesh, facet_domains, ind, c)
    #r = np.sqrt(A / pi) # This old estimate is a few % lower because of boundary discretization errors

    return A, c, r, normal

def compute_area(mesh, ind, facet_domains):
    # Some convenient variables
    assert facet_domains is not None
    dsi = ds[facet_domains](ind)

    # Compute area of boundary tesselation by integrating 1.0 over all facets
    A = assemble(Constant(1.0)*dsi, mesh=mesh)
    assert A > 0.0, "Expecting positive area, probably mismatch between mesh and markers!"
    return A

def compute_transient_scale_value(bc, period, mesh, facet_domains, ind, scale_value):
    dsi = ds[facet_domains](ind)
    form = sqrt(as_vector(bc)**2) * dsi

    N = 100
    qt = [0]*N
    for i, t in enumerate(np.linspace(0, period, N)):
        for e in bc:
            e.set_t(t)
        qt[i] = assemble(form, mesh=mesh)
    for e in bc:
        e.set_t(0.0)

    q_avg = sum(qt) / len(qt)
    return scale_value / q_avg