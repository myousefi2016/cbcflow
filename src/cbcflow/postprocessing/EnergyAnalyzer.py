
from .PPField import PPField
from math import sqrt
from dolfin import assemble

# FIXME: Split into separate fields
class EnergyAnalyzer(PPField):
    def compute(self, pp, spaces, problem):
        u = pp.get("Velocity")
        t = pp.get("t")

        dx = problem.dx

        u0_norm = sqrt(assemble(u[0]**2*dx()))
        u1_norm = sqrt(assemble(u[1]**2*dx()))
        u2_norm = sqrt(assemble(u[2]**2*dx()))
        energy = sqrt(u0_norm**2 +  u1_norm**2 + u2_norm**2)

        # FIXME: Dont print from fields, make it easy to control printing with parameters instead
        print "t ", t, " kinetic energy  ", energy

        data = {"u0" : u0_norm, "u1": u1_norm, "u2": u2_norm, "kinetic energy": energy}
        return data
