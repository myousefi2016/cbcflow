#!/usr/bin/env python

from cbcflow import *
from cbcflow.dol import *

from numpy import array

LENGTH = 4.0
RADIUS = 0.5

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > LENGTH*(1.0 - DOLFIN_EPS)

class Channel3D(NSProblem):
    "3D elongated box test problem with stationary flow state but no known analytical solution."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Create mesh
        N = self.params.N
        M = int(N*LENGTH/(2*RADIUS) + 0.5)
        mesh = UnitCubeMesh(M, N, N)
        x = mesh.coordinates()[:,0]
        y = mesh.coordinates()[:,1]
        z = mesh.coordinates()[:,2]
        x = LENGTH*x
        y = RADIUS*2*(y - 0.5)
        z = RADIUS*2*(z - 0.5)
        mesh.coordinates()[:,0] = x
        mesh.coordinates()[:,1] = y
        mesh.coordinates()[:,2] = z

        # We will apply markers with these id values
        self.wall_boundary_id = 0
        self.left_boundary_id = 1
        self.right_boundary_id = 2
        self.undefined_boundary_id = 3

        # Create boundary markers
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(self.undefined_boundary_id)
        DomainBoundary().mark(facet_domains, self.wall_boundary_id)
        Left().mark(facet_domains, self.left_boundary_id)
        Right().mark(facet_domains, self.right_boundary_id)

        # Setup analytical solution constants
        self.Upeak = 1.0
        self.U = self.Upeak / RADIUS**4
        nu = self.params.mu / self.params.rho
        self.beta = 2.0 * nu * self.U # FIXME: Can we estimate this for a square channel? This is for a pipe.

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=30e-5,
            dt=1e-5,
            # Physical parameters
            rho=1.0,
            mu=1.0/10.0,
            )
        params.update(
            # Spatial parameters
            N=6,
            )
        return params

    def analytical_solution(self, spaces, t):
        ux = Expression("U*(radius*radius - x[1]*x[1])*(radius*radius - x[2]*x[2])", U=1.0, radius=RADIUS)
        ux.U = self.U
        c0 = Constant(0)
        u0 = [ux, c0, c0]

        p0 = Expression("-beta*x[0]", beta=1.0)
        p0.beta = self.beta

        return (u0, p0)

    def initial_conditions(self, spaces, controls):
        return self.analytical_solution(spaces, 0.0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        ua, pa = self.analytical_solution(spaces, t)

        # Create no-slip and inflow boundary condition for velocity
        c0 = Constant(0)
        bcu = [
            ([c0, c0, c0], self.wall_boundary_id),
            (ua, self.left_boundary_id),
            ]

        # Create outflow boundary conditions for pressure
        bcp = [(pa, self.right_boundary_id)]

        return (bcu, bcp)

def main():
    problem = Channel3D()
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
