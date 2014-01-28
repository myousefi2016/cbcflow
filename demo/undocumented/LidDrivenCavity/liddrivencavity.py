#!/usr/bin/env python
__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-03-19"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"
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
            dt=0.01,
            T=2.5,
            # Physical parameters
            rho=1.0,
            mu=1.0/1000.0,
            )
        params.update(
            # Spatial parameters
            N=32,
            )
        return params
    
    def test_functionals(self):
        return [Minimum("StreamFunction", {"save": False, "start_time": 2.5-DOLFIN_EPS, "end_time": 2.5+DOLFIN_EPS})]
        
    def test_references(self):
        return [-0.061076605]
        

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


def main():
    set_log_level(100)
    problem = LidDrivenCavity({"N": 256})
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
    scheme = Yosida({"u_degree": 2})
    #scheme = IPCS()
    #scheme = IPCS({"u_degree": 2})
    #scheme = IPCS_Stable({"theta": 0.5})
    
    #solver_p_neumann=("gmres", "hypre_amg"),
    postprocessor = NSPostProcessor({"casedir": "Results_"+str(scheme.shortname())})
    
    #velocity = Velocity({"save": True})
    #pressure = Pressure({"save": True})
    #psi = StreamFunction({"save": True})
    #cfl = LocalCfl({"save": True})
    
    #postprocessor.add_fields([velocity, pressure, psi])
    test_functionals = problem.test_functionals()
    postprocessor.add_fields(test_functionals)
    
    solver = NSSolver(problem, scheme, postprocessor)
    ns = solver.solve()
    
    num_dofs = ns["spaces"].V.dim()+ns["spaces"].Q.dim()
    print num_dofs
    #print ns
    
    for i, tf in enumerate(test_functionals):
        val = postprocessor.get(tf.name)
        ref = problem.test_references()[i]
        err = sqrt((val-ref)**2)
        rel_err = err/abs(ref)
        
        print tf.name, val, ref, err, rel_err

if __name__ == "__main__":
    main()
