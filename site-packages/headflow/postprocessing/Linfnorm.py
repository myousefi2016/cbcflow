from .MetaPPField import MetaPPField
from dolfin import *

class Linfnorm(MetaPPField):
    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)
        if isinstance(u, float):
            return abs(u)
        else:
            return u.vector().norm("linf")
