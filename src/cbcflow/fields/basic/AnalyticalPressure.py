from ..bases.PPField import PPField
from dolfin import as_vector, Function

class AnalyticalPressure(PPField):
    @classmethod
    def default_params(cls):
        params = PPField.default_params()
        params.replace(
            assemble=False,
            project=True,
            interpolate=False,
            )
        # TODO: Perhaps we should require that analytical_solution returns an Expression or Function and use interpolate instead?
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

        return self.expr2function(pa, self._function)
