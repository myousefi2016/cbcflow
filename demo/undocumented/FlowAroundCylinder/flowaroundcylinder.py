#!/usr/bin/env python
__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-03-19"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from cbcflow import *
from cbcflow.dol import *

class LeftBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0)

class RightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 10.0)

class Cylinder(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (sqrt((x[0]-2.0)**2+(x[1]-0.5)**2) < 0.12+DOLFIN_EPS)

class FlowAroundCylinder(NSProblem):
    "Flow around a cylinder in 2D."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Create mesh
        r = Rectangle(0,0, 10, 1)
        c = Circle(2.0, 0.5, 0.12)
        mesh = Mesh(r-c, self.params.N)

        # Create boundary markers
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(4)
        DomainBoundary().mark(facet_domains, 0)
        Cylinder().mark(facet_domains, 0)
        LeftBoundary().mark(facet_domains, 1)
        RightBoundary().mark(facet_domains, 2)

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=5.0,
            dt=0.1,
            # Physical parameters
            rho=1.0,
            mu=1.0/1000.0,
            )
        params.update(
            # Spatial parameters
            N=64,
            )
        return params

    def initial_conditions(self, spaces, controls):
        c0 = Constant(0)
        u0 = [c0, c0]
        p0 = c0
        return (u0, p0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        c0 = Constant(0)
        c1 = Constant(1)

        # Create no-slip boundary condition for velocity
        bcu1 = ([c1, c0], 1)
        bcu2 = ([c0, c0], 0)
        bcu3 = ([c1, c0], 2)

        # Create boundary conditions for pressure
        bcp1 = (c0, 2)

        # Collect and return
        bcu = [bcu1, bcu2]
        bcp = []
        if 1:
            bcp += [bcp1]
        else:
            bcu += [bcu3]
        return (bcu, bcp)

    '''
    # FIXME: Change this to use the new test_functionals, test_references interface:
    def functional(self, t, u, p):
        # Only check final time
        if t < self.T:
            return 0
        else:
            # Compute stream function and report minimum
            psi = StreamFunction(u)
            vals  = psi.vector().array()
            vmin = MPI.min(vals.min())

            cbcflow_print("Stream function has minimal value %s" % vmin)

            return vmin

    def reference(self, t):
        # Only check final time
        if t < self.T:
            return 0.0
        return -0.061076605
    '''

def main():
    problem = FlowAroundCylinder()
    scheme = IPCS_Stable()

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
