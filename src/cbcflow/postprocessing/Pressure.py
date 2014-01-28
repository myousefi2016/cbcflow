from .PPField import PPField
from dolfin import Function, FunctionAssigner
from ..core.spaces import NSSpacePoolMixed, NSSpacePoolSplit, NSSpacePoolSegregated

class Pressure(PPField):
    def convert(self, pp, spaces, problem):
        # Hack to get given p in whatever format it has
        p = super(Pressure, self).convert(pp, spaces, problem)

        if not isinstance(p, Function):
            if not hasattr(self, "_p"):
                self._p = Function(spaces.Q)
                assert isinstance(spaces, NSSpacePoolMixed)
                self._assigner = FunctionAssigner(spaces.Q, spaces.W.sub(1))

            # Hack: p is a Indexed(Coefficient()),
            # get the underlying mixed function
            w = p.operands()[0]
            self._assigner.assign(self._p, w.sub(1))

            p = self._p

        assert isinstance(p, Function)
        return p

SolverPressure = Pressure
