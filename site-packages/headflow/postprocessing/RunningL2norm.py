from .MetaPPField import MetaPPField
from dolfin import *
import numpy as np

class RunningL2norm(MetaPPField):
    def before_first_compute(self, pp, spaces, problem):
        self._count = 0
        self._sum = 0
        self._value = 0

    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)

        self._count += 1
        self._sum += u**2
        self._value = np.sqrt(self._sum) / self._count

        return self._value
