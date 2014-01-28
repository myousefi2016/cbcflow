from ..bases.MetaPPField import MetaPPField
from dolfin import Function, MPI
import numpy

class Maximum(MetaPPField):
    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)
        
        if isinstance(u, Function):
            return MPI.max(numpy.max(u.vector().array()))
        else:
            return MPI.max(max(u))
        