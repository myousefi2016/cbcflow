from ..bases.PPField import PPField
from dolfin import *

class VelocityDivergence(PPField):
    def compute(self, pp, spaces, problem):
        u = pp.get("Velocity")
        return sqrt(assemble(div(u)**2*dx()))
