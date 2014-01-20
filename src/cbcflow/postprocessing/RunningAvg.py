from .MetaPPField import MetaPPField
from dolfin import *
import numpy

class RunningAvg(MetaPPField):
    def before_first_compute(self, pp, spaces, problem):
        self._value = None
        self._count = 0

    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)

        self._count += 1

        if self._value is None:
            if isinstance(u, Function):
                self._sum = Function(u)
                self._value = Function(u)
            else:
                self._sum = u
                self._value = u
        else:
            if isinstance(u, Function):
                self._sum.vector().axpy(1.0, u.vector())
                self._value.vector().zero()
                self._value.vector().axpy(1.0/self._count, self._sum.vector())
            else:
                self._sum += u
                self._value = u / self._count

        return self._value

    def after_last_compute(self, pp, spaces, problem):
        return self._value
