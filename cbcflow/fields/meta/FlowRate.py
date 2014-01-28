
from ..bases.PPField import PPField
from dolfin import *

class FlowRate(PPField):
    def __init__(self, boundary_id, params=None, label=None):
        PPField.__init__(self, params, label)
        self.boundary_id = boundary_id

    @property
    def name(self):
        n = "%s_%s" % (self.__class__.__name__, self.boundary_id)
        if self.label: n += "_"+self.label
        return n

    def compute(self, pp, spaces, problem):
        u = pp.get("Velocity")

        n = problem.mesh.ufl_cell().n
        dsi = problem.ds(self.boundary_id)

        M = dot(u, n)*dsi
        return assemble(M)
