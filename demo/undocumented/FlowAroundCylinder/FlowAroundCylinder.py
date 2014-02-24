#!/usr/bin/env python

from cbcflow import *
from cbcflow.dol import *

from os import path

files = [path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/cylinder_0.6k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/cylinder_2k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/cylinder_8k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/cylinder_32k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/cylinder_129k.xml.gz"),
        ]

class LeftBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0)

class RightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 10.0)

class Cylinder(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (sqrt((x[0]-2.0)**2+(x[1]-0.5)**2) < 0.12+DOLFIN_EPS)

class Wall(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[1], 0.0) or near(x[1], 1.0))

class FlowAroundCylinder(NSProblem):
    "Flow around a cylinder in 2D."
    
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
            refinement_level=0,
            )
        return params

    def __init__(self, params=None):
        NSProblem.__init__(self, params)
        
        # Load mesh
        mesh = Mesh(files[self.params.refinement_level])

        # Create boundary markers
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(4)
        Wall().mark(facet_domains, 0)
        Cylinder().mark(facet_domains, 0)
        LeftBoundary().mark(facet_domains, 1)
        RightBoundary().mark(facet_domains, 2)

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

    def initial_conditions(self, spaces, controls):
        c0 = Constant(0)
        u0 = [c0, c0]
        p0 = c0
        return (u0, p0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        c0 = Constant(0)
        c1 = Constant(1)

        # Create no-slip boundary condition for velocity
        bcu0 = ([c0, c0], 0)
        bcu1 = ([c1, c0], 1)

        # Create boundary conditions for pressure
        bcp0 = (c0, 2)

        # Collect and return
        bcu = [bcu0, bcu1]
        bcp = [bcp0]
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
    problem = FlowAroundCylinder({"refinement_level": 2})
    scheme = IPCS_Stable()

    casedir = "results_demo_%s_%s" % (problem.shortname(), scheme.shortname())
    postprocessor = NSPostProcessor({"casedir": casedir})
    
    plot_and_save = dict(plot=True, save=True)
    fields = [
        Pressure(plot_and_save),
        Velocity(plot_and_save),
        StreamFunction(plot_and_save),
        ]
    postprocessor.add_fields(fields)

    solver = NSSolver(problem, scheme, postprocessor)
    solver.solve()

if __name__ == "__main__":
    main()
