from .PPField import PPField
from dolfin import *

class VelocityCurl(PPField):
    def compute(self, pp, spaces, problem):
        u = pp.get("Velocity")
        return sqrt(assemble(curl(u)**2*dx()))
