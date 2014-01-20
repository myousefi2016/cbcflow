from .PPField import PPField
from ..core.utils import sigma
from dolfin import *

class Stress(PPField):
    def before_first_compute(self, pp, spaces, problem):
        if self.params.assemble:
            V = spaces.get_space(0, 2)
        else:
            V = spaces.DV
        self._function = Function(V, name=self.name)

    def compute(self, pp, spaces, problem):
        u = pp.get("Velocity")
        p = pp.get("Pressure")
        mu = problem.params.mu

        expr = sigma(u, p, mu)
        #u*epsilon(u) - p*Identity(u.cell().d) # TODO: is this with negative pressure?

        return self.expr2function(expr, self._function)
