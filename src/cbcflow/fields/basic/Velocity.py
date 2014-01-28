from ..bases.PPField import PPField
from dolfin import Function, FunctionAssigner, error
from ...core.spaces import NSSpacePoolMixed, NSSpacePoolSplit, NSSpacePoolSegregated

class Velocity(PPField):
    def convert(self, pp, spaces, problem):
        # Hack to get given u in whatever format it has,
        # avoiding circular reference to this field
        u = super(Velocity, self).convert(pp, spaces, problem)
        d = spaces.d

        if not isinstance(u, Function):
            if not hasattr(self, "_u"):
                self._u = Function(spaces.V)

                if isinstance(spaces, NSSpacePoolMixed):
                    self._assigner = FunctionAssigner(spaces.V, spaces.W.sub(0))
                elif isinstance(spaces, NSSpacePoolSegregated):
                    self._assigner = FunctionAssigner(spaces.V, [spaces.U]*d)
                else:
                    error("It doesnt make sense to create a function assigner for a split space.")

            if isinstance(spaces, NSSpacePoolMixed):
                # Hack: u is a ListTensor([Indexed(Coefficient()),...]),
                # get the underlying mixed function
                w = u.operands()[0].operands()[0]
                assert w.shape() == (d+1,)
                us = w.sub(0)

            elif isinstance(spaces, NSSpacePoolSegregated):
                us = [u[i] for i in range(d)]

            self._assigner.assign(self._u, us)
            u = self._u

        assert isinstance(u, Function)
        return u
