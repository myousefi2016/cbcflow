from .MetaPPField import MetaPPField
from dolfin import *
import numpy

class RunningMin(MetaPPField):
    def before_first_compute(self, pp, spaces, problem):
        self._value = None

    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)

        if self._value is None:
            if isinstance(u, Function):
                self._value = Function(u)
            else:
                self._value = u
        else:
            if isinstance(u, Function):
                # TODO: Test! This might work at least in serial, what about paralell?
                self._value.vector()[:] = numpy.min(self._value.vector()[:], u.vector()[:])
            else:
                self._value = min(self._value, u)

        return self._value

    def after_last_compute(self, pp, spaces, problem):
        return self._value
