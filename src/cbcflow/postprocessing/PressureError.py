from .PPField import PPField
from dolfin import as_vector, Function

class PressureError(PPField):
    @classmethod
    def default_params(cls):
        params = PPField.default_params()
        params.replace(
            assemble=False,
            project=True,
            interpolate=False,
            )
        return params

    def before_first_compute(self, pp, spaces, problem):
        if self.params.assemble:
            degree = 0
        else:
            degree = spaces.p_degree + 1 # TODO: Is +1 sufficient?
        V = spaces.get_space(degree, 0)
        self._function = Function(V, name=self.name)

    def compute(self, pp, spaces, problem):
        p = pp.get("Pressure")
        t = pp.get("t")

        ua, pa = problem.analytical_solution(spaces, t)
        pe = pa - p

        return self.expr2function(pe, self._function)
