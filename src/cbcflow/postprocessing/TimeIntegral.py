from .MetaPPField import MetaPPField
from dolfin import *

class TimeIntegral(MetaPPField):
    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)

        t1 = pp.get("t")
        t0 = pp.get("t", -1)
        dt = t1 - t0

        if isinstance(u, Function):
            # Create placeholder for sum the first time
            if not hasattr(self, "_sum"):
                self._sum = Function(u.function_space())
            # Accumulate using backward euler integration
            self._sum.vector().axpy(dt, u.vector()) # FIXME: Validate this, not tested!
        else:
            # Create placeholder for sum the first time
            if not hasattr(self, "_sum"):
                self._sum = 0.0
            # Accumulate using backward euler integration
            self._sum += u

        # Increase count and return sum
        if not hasattr(self, "_count"):
            self._count = 0
        self._count += 1
        return self._sum
