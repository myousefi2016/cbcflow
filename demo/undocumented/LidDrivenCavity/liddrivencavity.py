#!/usr/bin/env python
__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-03-19"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"
import sys
sys.path.insert(0, '../../../src/')
from cbcflow import *
from dolfin import *

set_log_level(100)

class Lid(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] > DOLFIN_EPS and x[0] < 1.0 - DOLFIN_EPS and x[1] > 1.0 - DOLFIN_EPS


class LidDrivenCavity(NSProblem):
    "2D lid-driven cavity test problem."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)
        N = self.params.N
        mesh = UnitSquareMesh(N, N)

        # Create boundary markers
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(2)
        DomainBoundary().mark(facet_domains, 0)
        Lid().mark(facet_domains, 1)

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            dt=0.001,
            T=0.5,
            # Physical parameters
            rho=1.0,
            mu=1.0/10.0,
            )
        params.update(
            # Spatial parameters
            N=32,
            )
        return params

    def initial_conditions(self, spaces, controls):
        u0 = [Constant(0), Constant(0)]
        p0 = Constant(0)
        return (u0, p0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        # Set "lid" in motion, no-slip condition elsewhere
        bcu = [([Constant(0), Constant(0)], 0),
               ([Constant(1), Constant(0)], 1)]
        bcp = []
        return (bcu, bcp)


if __name__ == "__main__":
    #set_log_level(0)
    problem = LidDrivenCavity({"N": 16})
    #parameters["krylov_solver"]["relative_tolerance"] = 1e-15
    parameters["linear_algebra_backend"] = "PETSc"
    #print parameters["krylov_solver"].items()
    #exit()
    
    #show_problem(problem)
    #scheme = IPCS({"solver_p_neumann": ("cg",)})
    #print dir(scheme)
    #print scheme.shortname()
    #exit()
    #scheme = IPCS({"theta": 0.5, "solver_p_neumann": ("gmres", "jacobi")})
    #scheme = IPCS({"theta": 0.5})
    scheme = Yosida()
    #scheme = IPCS()
    #scheme = IPCS({"u_degree": 2})
    #scheme = IPCS_Stable({"theta": 1.0})
    
    #solver_p_neumann=("gmres", "hypre_amg"),
    postprocessor = NSPostProcessor({"casedir": "Results_"+str(scheme.shortname())})
    
    velocity = Velocity({"save": True, "save_as": 'pvd'})
    pressure = Pressure({"save": True, "save_as": 'pvd'})
    #cfl = LocalCfl({"save": True})
    
    postprocessor.add_fields([velocity, pressure])
    
    solver = NSSolver(problem, scheme, postprocessor)
    solver.solve()
    
    
    
    #from demo_main import demo_main
    #demo_main(LidDrivenCavity)
