from ..dol import *

def compute_resistance_value(C, u, ind, facet_domains):
    mesh = facet_domains.mesh()
    dsi = ds[facet_domains](ind)
    n = FacetNormal(mesh)
    u = as_vector(u)
    form = dot(u, n)*dsi
    return C*assemble(form)

# TODO: [martin] don't understand this design, why subclassing Constant?
class Resistance(Constant):
    def __init__(self, C, u, ind, facet_domains):
        Constant.__init__(self, compute_resistance_value(C, u, ind, facet_domains))
