from ..bases.PPField import PPField
from dolfin import assemble, dx, sqrt, TestFunction, Function, Constant

class LocalCfl(PPField):
    def before_first_compute(self, pp, spaces, problem):
        DG0 = spaces.get_space(0, 0)
        self._v = TestFunction(DG0)
        self._cfl = Function(DG0)

    def compute(self, pp, spaces, problem):
        t1 = pp.get("t")
        t0 = pp.get("t", -1)
        dt = Constant(t1 - t0)
        u = pp.get("Velocity")

        cell = problem.mesh.ufl_cell()
        hF = cell.circumradius
        hK = cell.volume
        scaling = 1.0 / hK

        assemble((dt * sqrt(u**2) / hF)*self._v*scaling*dx(),
                 tensor=self._cfl.vector())

        return self._cfl
