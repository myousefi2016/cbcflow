
from ..core.paramdict import ParamDict
from .PPField import PPField

from dolfin import FunctionSpace, TrialFunction, TestFunction, Function, grad, det, Constant

class Q(PPField):
    @classmethod
    def default_params(cls):
        params = PPField.default_params()
        params.replace(
            assemble=False, # Change to this to use assemble into DG0 by default
            project=True,
            interpolate=False,
            )
        return params

    def before_first_compute(self, pp, spaces, problem):
        if self.params.assemble:
            V = spaces.get_space(0, 0)
        else:
            # Accurate degree is 2*(spaces.u_degree-1)
            degree = 1
            V = spaces.get_space(degree, 0)
        self._function = Function(V, name=self.name)

    def compute(self, pp, spaces, problem):
        u = pp.get("Velocity")

        S = (grad(u) + grad(u).T)/2
        Omega = (grad(u) - grad(u).T)/2
        expr = Constant(0.5) * (Omega**2 - S**2)

        return self.expr2function(expr, self._function)
