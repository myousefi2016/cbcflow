from .PPField import PPField
from math import sqrt
from dolfin import assemble

class KineticEnergy(PPField):
    def compute(self, pp, spaces, problem):
        u = pp.get("Velocity")
        dx = problem.dx
        u_norms = [assemble(u[d]**2*dx()) for d in range(u.shape()[0])]
        energy = sqrt(sum(u_norms[d] for d in range(u.shape()[0])))
        return energy
