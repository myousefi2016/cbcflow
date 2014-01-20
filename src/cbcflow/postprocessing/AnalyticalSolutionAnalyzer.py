
from .PPField import PPField
from math import sqrt
from dolfin import assemble, errornorm, Constant

class AnalyticalSolutionAnalyzer(PPField):
    def compute(self, pp, spaces, problem):
        p = pp.get("Pressure")
        u = pp.get("Velocity")
        t = pp.get("t")

        mesh = problem.mesh
        dx = problem.dx

        u_anal, p_anal = problem.analytical_solution(spaces, t)

        u0_error = sqrt(assemble((u_anal[0] - u[0])**2*dx()))
        u1_error = sqrt(assemble((u_anal[1] - u[1])**2*dx()))
        u2_error = sqrt(assemble((u_anal[2] - u[2])**2*dx()))
        c = assemble((p_anal-p)*dx()) / assemble(Constant(1)*dx(), mesh=mesh)

        p.vector()[:] += c # FIXME FIXME FIXME Kent, what is this? You cannot modify the input functions here...
        p_error = errornorm(p_anal, p)
        
        p_dx = assemble(p*dx())
        p_anal_dx = assemble(p_anal*dx(), mesh=mesh)

        data = {
            #"u": u_error, # L2norm("VelocityError", ...)
            "u0" : u0_error,
            "u1": u1_error,
            "u2": u2_error,
            "p": p_error, # L2norm("PressureError", ...)
            "pdx": p_dx, # DomainAverage("Pressure", ...)
            "analytical_pdx" : p_anal_dx, # DomainAverage("AnalyticalPressure", ...)
            }
        #self.set_data(t, timestep, data)


        # FIXME: Dont print from fields, make it easy to control printing with parameters instead
        print "t ", t, " error (u0, u1, u2, p, pdx, analytical_pdx) ", u0_error, u1_error, u2_error, p_error, p_dx, p_anal_dx

        return data
