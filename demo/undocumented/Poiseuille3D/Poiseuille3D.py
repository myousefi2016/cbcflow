#!/usr/bin/env python

from cbcflow import *
from cbcflow.dol import *
from os import path

import numpy as np

files = [path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/pipe_1k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/pipe_3k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/pipe_24k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/pipe_203k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/pipe_1611k.xml.gz"),
        ]

class Inflow(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < 1e-6 and on_boundary

class Outflow(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 10.0-1e-6 and on_boundary


LENGTH = 10.0
RADIUS = 0.5

class Poiseuille3D(NSProblem):
    "3D pipe test problem with known stationary analytical solution."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Load mesh
        mesh = Mesh(files[self.params.refinement_level])

        # We know that the mesh contains markers with these id values
        self.wall_boundary_id = 0
        self.left_boundary_id = 1
        self.right_boundary_id = 2
        
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(3)
        DomainBoundary().mark(facet_domains, self.wall_boundary_id)
        Inflow().mark(facet_domains, self.left_boundary_id)
        Outflow().mark(facet_domains, self.right_boundary_id)
        
        # Setup analytical solution constants
        Q = self.params.Q

        #print "Using stationary bcs. Analytical solution should hold."
        self.Q_coeffs = [(0.0, Q), (1.0, Q)]

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=None,
            dt=1e-3,
            period=0.8,
            num_periods=0.1,
            # Physical parameters
            rho=10.0,
            mu=1.0/30.0,
            )
        params.update(
            # Spatial parameters
            refinement_level=0,
            # Analytical solution parameters
            Q=1.0,
            )
        return params

    def analytical_solution(self, spaces, t):
        A = pi*RADIUS**2
        Q = self.params.Q
        mu = self.params.mu       
        dpdx = 8.0*Q*mu/(A*RADIUS**2)
        
        ux = Expression("(radius*radius - x[1]*x[1] - x[2]*x[2])/(4*mu)*dpdx", radius=RADIUS, mu=mu, dpdx=dpdx, cell=tetrahedron, degree=2)
        uy = Constant(0.0)
        
        uz = Constant(0.0)
        u = [ux, uy, uz]
        
        p = Expression("dpdx * (length-x[0])", dpdx=dpdx, length=LENGTH)
        
        return (u,p)

    def test_references(self, spaces, t):
        return self.analytical_solution(spaces, t)
    
    def test_fields(self):
        return [SolutionField("Velocity"), SolutionField("Pressure")]

    def initial_conditions(self, spaces, controls):
        return self.analytical_solution(spaces, 0.0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        # Create no-slip bcs
        d = len(u)
        u0 = [Constant(0.0)] * d
        noslip = (u0, self.wall_boundary_id)

        # Create Poiseuille inflow bcs
        uin = make_poiseuille_bcs(self.Q_coeffs, self.mesh, self.left_boundary_id, self.params.Q, self.facet_domains)
        for ucomp in uin:
            ucomp.set_t(t)
        inflow = (uin, self.left_boundary_id)

        # Create outflow bcs for pressure
        _, pa = self.analytical_solution(spaces, t)
        outflow = (pa, self.right_boundary_id)

        # Return bcs in two lists
        bcu = [noslip, inflow]
        bcp = [outflow]
        return (bcu, bcp)

def main():
    problem = Poiseuille3D({"refinement_level": 0})
    scheme = IPCS_Stable({"u_degree": 1})

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
