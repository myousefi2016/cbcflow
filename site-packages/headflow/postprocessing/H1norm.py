from .MetaPPField import MetaPPField
from dolfin import *

class H1norm(MetaPPField):
    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)
        M = u**2*dx() + grad(u)**2*dx()
        return sqrt(assemble(M))
