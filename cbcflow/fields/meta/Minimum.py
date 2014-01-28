from ..bases.MetaPPField import MetaPPField
from dolfin import Function, MPI
import numpy

class Minimum(MetaPPField):
    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)
        
        if isinstance(u, Function):
            return MPI.min(numpy.min(u.vector().array()))
        else:
            return MPI.min(min(u))
        