
from ..bases.MetaPPField2 import MetaPPField2
from dolfin import *

class DiffH1norm(MetaPPField2):
    "Compute the full H1 norm of the difference between uh and u relative to u."
    def compute(self, pp, spaces, problem):
        uh = pp.get(self.valuename1)
        u = pp.get(self.valuename2)
        e = uh - u

        dx = problem.dx

        uh2 = assemble((uh**2 + grad(uh)**2)*dx(), mesh=problem.mesh)
        u2 = assemble((u**2 + grad(u)**2)*dx(), mesh=problem.mesh)
        e2 = assemble((e**2 + grad(e)**2)*dx(), mesh=problem.mesh)

        # Relative norm
        eps = 1e-14
        if abs(u2) > eps:
            return sqrt(e2 / u2)
        else:
            return sqrt(e2)
