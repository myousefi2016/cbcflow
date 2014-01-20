from .MetaPPField import MetaPPField, MetaPPField2
from dolfin import *

class DiffH1seminorm(MetaPPField2):
    "Compute the H1 semi norm of the difference between uh and u relative to u."
    def compute(self, pp, spaces, problem):
        uh = pp.get(self.valuename1)
        u = pp.get(self.valuename2)
        e = uh - u

        dx = problem.dx

        uh2 = assemble((grad(uh)**2)*dx(), mesh=problem.mesh)
        u2 = assemble((grad(u)**2)*dx(), mesh=problem.mesh)
        e2 = assemble((grad(e)**2)*dx(), mesh=problem.mesh)

        # Relative norm
        eps = 1e-14
        if abs(u2) > eps:
            return sqrt(e2 / u2)
        else:
            return sqrt(e2)
