from ..dol import *

from .bcutils import compute_area

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
