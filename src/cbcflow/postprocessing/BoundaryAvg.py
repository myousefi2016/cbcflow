from .MetaPPField import MetaPPField
from dolfin import assemble, ds

class BoundaryAvg(MetaPPField):
    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)
        value = assemble(u*ds(), mesh=problem.mesh)
        return value
