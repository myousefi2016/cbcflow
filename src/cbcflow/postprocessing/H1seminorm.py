from .MetaPPField import MetaPPField
from dolfin import *

class H1seminorm(MetaPPField):
    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)
        M = grad(u)**2*dx()
        return sqrt(assemble(M))
