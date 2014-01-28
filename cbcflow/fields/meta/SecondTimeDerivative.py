from ..bases.MetaPPField import MetaPPField
from dolfin import *

class SecondTimeDerivative(MetaPPField):
    def compute(self, pp, spaces, problem):
        u2 = pp.get(self.valuename)
        u1 = pp.get(self.valuename, -1)
        u0 = pp.get(self.valuename, -2)

        t2 = pp.get("t")
        t1 = pp.get("t", -1)
        t0 = pp.get("t", -2)
        dt1 = t2 - t1
        dt0 = t1 - t0

        # Computing d2u/dt2 as if t1 = 0.5*(t2+t0), i.e. assuming fixed timesteps
        # TODO: Find a more accurate formula if dt1 != dt0?
        dt = 0.5*(dt1 + dt0)

        if isinstance(u0, Function):
            # Create function to hold result first time,
            # assuming u2, u1 and u0 are Functions in same space
            if not hasattr(self, "_du"):
                self._du = Function(u0.function_space())

            # Compute finite difference derivative # FIXME: Validate this, not tested!
            self._du.vector().zero()
            self._du.vector().axpy(+1.0/dt**2, u2.vector())
            self._du.vector().axpy(-2.0/dt**2, u1.vector())
            self._du.vector().axpy(+1.0/dt**2, u0.vector())
            return self._du
        else:
            # Assuming scalar value
            du = (u2 - 2.0*u1 + u0) / dt**2
            return du
