from .PPField import PPField
from dolfin import Function
from ..core.spaces import NSSpacePoolMixed, NSSpacePoolSplit, NSSpacePoolSegregated

class PhysicalPressure(PPField):
    "The physical pressure is the solver pressure scaled by density."
    def compute(self, pp, spaces, problem):
        # Get solver pressure
        p = pp.get("Pressure")

        if not hasattr(self, "_p"):
            self._p = Function(spaces.Q)

        if isinstance(p, Function):
            self._p.assign(p)
        else:
            assert isinstance(spaces, NSSpacePoolMixed)
            if not hasattr(self, "_assigner"):
                self._assigner = FunctionAssigner(spaces.Q, spaces.W.sub(1))
            # Hack: p is a Indexed(Coefficient()),
            # get the underlying mixed function
            w = p.operands()[0]
            self._assigner.assign(self._p, w.sub(1))

        # Scale by density
        rho = problem.params.rho
        pv = self._p.vector()
        pv *= rho
        p = self._p

        assert isinstance(p, Function)
        return p
