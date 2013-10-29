from .MetaPPField import MetaPPField
from dolfin import *
import numpy

class Magnitude(MetaPPField):
    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)
        
        if isinstance(u, Function):
            if u.rank() == 0:
                return u
            elif u.rank() >= 1:
                # Assume all subpaces are equal
                V = u.function_space().extract_sub_space([0]).collapse()
                mag = project(sqrt(inner(u,u)), V)
                return mag
        else:
            # Don't know how to handle object
            headlow_warning("Don't know how to calculate magnitude of object of type %s. Returning object." %type(u))
            return u
                        
    