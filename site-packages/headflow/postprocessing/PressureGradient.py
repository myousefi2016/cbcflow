from .PPField import PPField
from dolfin import *

class PressureGradient(PPField):
    def before_first_compute(self, pp, spaces, problem):
        if self.params.assemble:
            V = spaces.get_space(0, 1)
        else:
            V = spaces.DQ
        self._function = Function(V, name=self.name)

    def compute(self, pp, spaces, problem):
        p = pp.get("Pressure")

        expr = grad(p)

        return self.expr2function(expr, self._function)
