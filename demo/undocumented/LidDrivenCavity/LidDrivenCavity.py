#!/usr/bin/env python
from cbcflow import *
from dolfin import *

class Lid(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] > DOLFIN_EPS and x[0] < 1.0 - DOLFIN_EPS and x[1] > 1.0 - DOLFIN_EPS

class LidDrivenCavity(NSProblem):
    "2D lid-driven cavity test problem."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)
        
        refinements = [16,32,64,128,256]
        
        #N = self.params.N
        N = refinements[self.params.refinement_level]
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
            #N=32,
            refinement_level=0,
            )
        return params
    
    def test_fields(self):
        return [Minimum("StreamFunction", {"save": False, "start_time": 2.5-DOLFIN_EPS, "end_time": 2.5+DOLFIN_EPS})]
        
    def test_references(self, spaces, t):
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
    problem = LidDrivenCavity({"refinement_level": 1})
    
    
    from imp import find_module
    try:
        find_module("block")
        has_cbcblock = True
    except:
        has_cbcblock = False
    
    if has_cbcblock:
        scheme = Yosida() # Requires cbc.block
    else:
        scheme = IPCS_Stable({"solver_p_neumann": ("cg", "ilu")}) # Displays pressure oscillations

    casedir = "results_demo_%s_%s" % (problem.shortname(), scheme.shortname())
    plot_and_save = dict(plot=True, save=True)
    fields = [
        Pressure(plot_and_save),
        Velocity(plot_and_save),
        ]
    postproc = NSPostProcessor({"casedir": casedir})
    postproc.add_fields(fields)

    solver = NSSolver(problem, scheme, postproc)
    solver.solve()

if __name__ == "__main__":
    main()
