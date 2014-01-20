from .PPField import PPField
from dolfin import *

class Strain(PPField):
    def before_first_compute(self, pp, spaces, problem):
        if self.params.assemble:
            V = spaces.get_space(0, 2)
        else:
            V = spaces.DV
        self._function = Function(V, name=self.name)

    def compute(self, pp, spaces, problem):
        u = pp.get("Velocity")
        Du = grad(u)

        expr = 0.5*(Du + Du.T)

        return self.expr2function(expr, self._function)
