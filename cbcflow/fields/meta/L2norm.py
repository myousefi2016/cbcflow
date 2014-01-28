from ..bases.MetaPPField import MetaPPField
from dolfin import *

class L2norm(MetaPPField):
    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)
        if isinstance(u, float):
            return sqrt(u**2)
        else:
            M = u**2*dx()
            return sqrt(assemble(M))
