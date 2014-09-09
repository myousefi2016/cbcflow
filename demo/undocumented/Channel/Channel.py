#!/usr/bin/env python


from cbcflow import *
from cbcflow.dol import *

from numpy import array

class InflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS

class OutflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 1 - DOLFIN_EPS

class Channel(NSProblem):
    "2D channel test problem with known analytical solution."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Create mesh
        N = self.params.N
        mesh = UnitSquareMesh(N, N)

        # Create boundary markers
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(4)
        DomainBoundary().mark(facet_domains, 0)
        InflowBoundary().mark(facet_domains, 1)
        OutflowBoundary().mark(facet_domains, 2)

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=1.0e2/8,#0.5,
            dt=1.0/80,
            # Physical parameters
            rho=1.0,
            mu=1.0/8.0,
            )
        params.update(
            # Spatial parameters
            N=16,
            )
        return params

    def initial_conditions(self, spaces, controls):
        c0 = Constant(0)
        u0 = [c0, c0]
        p0 = Expression("1 - x[0]")
        return (u0, p0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        # Create no-slip boundary condition for velocity
        c0 = Constant(0)
        bcu = [([c0, c0], 0)]

        # Create boundary conditions for pressure
        p = Expression("1 - x[0]")
        bcp = [(p, 1),
               (p, 2)]

        return (bcu, bcp)

    """
    # FIXME: Change this to use the new test_functionals, test_references interface:
    def functional(self, t, u, p):
        if t < self.T:
            return 0
        else:
            return self.uEval(u, 0, (1.0, 0.5))

    def reference(self, t):
        if t < self.T:
            return 0
        else:
            num_terms = 10000
            u = 1.0
            c = 1.0
            for n in range(1, 2*num_terms, 2):
                a = 32.0 / (DOLFIN_PI**3*n**3)
                b = (1/8.0)*DOLFIN_PI**2*n**2
                c = -c
                u += a*exp(-b*t)*c
            return u

    def tolerance(self, problem):
        return 1e-11
    """

def main():
    problem = Channel()
    scheme = IPCS_Stable()

    casedir = "results_demo_%s_%s" % (problem.shortname(), scheme.shortname())
    plot_and_save = dict(plot=True, save=True)
    fields = [
        SolutionField("Pressure", plot_and_save),
        SolutionField("Velocity", plot_and_save),
        ]
    postproc = PostProcessor({"casedir": casedir})
    postproc.add_fields(fields)

    solver = NSSolver(problem, scheme, postproc)
    solver.solve()

if __name__ == "__main__":
    main()
