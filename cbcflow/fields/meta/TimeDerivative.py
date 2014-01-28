from ..bases.MetaPPField import MetaPPField
from dolfin import *

class TimeDerivative(MetaPPField):
    def compute(self, pp, spaces, problem):
        u1 = pp.get(self.valuename)
        u0 = pp.get(self.valuename, -1)

        t1 = pp.get("t")
        t0 = pp.get("t", -1)
        dt = t1 - t0

        if isinstance(u0, Function):
            # Create function to hold result first time,
            # assuming u1 and u0 are Functions in same space
            if not hasattr(self, "_du"):
                self._du = Function(u0.function_space())

            # Compute finite difference derivative # FIXME: Validate this, not tested!
            self._du.vector().zero()
            self._du.vector().axpy(+1.0/dt, u1.vector())
            self._du.vector().axpy(-1.0/dt, u0.vector())
            return self._du
        else:
            # Assuming scalar value
            du = (u1 - u0) / dt
            return du
