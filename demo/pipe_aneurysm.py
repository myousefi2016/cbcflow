#!/usr/bin/env python
__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2008-11-21"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Harish Narayanan, 2009.
# Modified by Anders Logg, 2010.
# Modified by Martin Alnaes, 2013.

#from scipy import *
from headflow import *
from headflow.dol import *

# Inflow boundary
class Inflow(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] < DOLFIN_EPS

# Outflow boundary
class Outflow(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] > 0.05 - DOLFIN_EPS

# Define the aneurysm region, everything outside the cylinder
class AneurysmCutoff(Expression):
    def eval(self, values, x):
        r = sqrt(x[1]**2 + x[2]**2)
        # FIXME: Is this well-defined?
        if r < 0.002 + 0.0001:
            values[0] = 0.0
        else:
            values[0] = 1.0

c0 = Constant(0)

class PipeAneurysm(NSProblem):
    "3D artery with a saccular aneurysm."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Load mesh
        if self.params.refinement_level > 4:
            raise RuntimeError("No mesh available for refinement level %d" % self.params.refinement_level)
        meshfilename = "../data/aneurysm_%d.xml.gz" % self.params.refinement_level
        mesh = Mesh(meshfilename)

        # Mark domains, 0 = noslip, 1 = inflow, 2 = outflow, 3 = rest
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(3)
        DomainBoundary().mark(facet_domains, 0)
        Inflow().mark(facet_domains, 1)
        Outflow().mark(facet_domains, 2)

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

    @classmethod
    def default_user_params(cls):
        params = ParamDict(
            # Spatial discretization parameters
            refinement_level = 1,
            # Time parameters
            dt = 0.005,
            T = 0.05,
            # Physical parameters
            rho = 1.0,
            mu = 3.5 / 1.025e6,
            )
        return params

    def initial_conditions(self, spaces, controls):
        u0 = [c0, c0, c0]
        p0 = c0
        return (u0, p0)

    def boundary_conditions(self, spaces, u, p, t, controls):

        # Create no-slip boundary condition for velocity
        g_noslip = [c0, c0, c0]
        bc_noslip = (g_noslip, 0)

        # Create inflow boundary condition for velocity
        inflow_exprs = ('1.0*(sin(30*t))*(1.0-(x[1]*x[1]+x[2]*x[2])/(r*r))',
                        '0.0*(sin(30*t))*(1.0-(x[1]*x[1]+x[2]*x[2])/(r*r))',
                        '0.0*(sin(30*t))*(1.0-(x[1]*x[1]+x[2]*x[2])/(r*r))')
        g_inflow = [Expression(e, degree=3, r=0.002, t=0.0) for e in inflow_exprs]
        for g in g_inflow: g.t = float(t)
        bc_inflow = (g_inflow, 1)

        # Create outflow boundary condition for pressure
        g_outflow = c0
        bc_outflow = (g_outflow, 2)

        # Collect boundary conditions
        bcu = [bc_noslip, bc_inflow]
        bcp = [bc_outflow]

        return (bcu, bcp)

    def update(self, spaces, u, p, t, timestep, bcs, observations, controls):
        bcu, bcp = bcs
        inflow = bcu[1][0]
        for ue in inflow: ue.t = float(t)

# Old code:
'''
    def functional(self, t, u, p):
         if t < self.T:
             return 0.0
         x = (0.025, -0.006, 0.0)
         # TODO: Is this the same as the original
         #       uEval(u, 0, x)?
         return parallell_eval(u[0], x)

    def reference(self, t):
        """The reference value was computed using on a fine mesh
        (level 6). Values obtained for different refinement levels
        are listed below for Chorin and IPCS.

              Chorin                 IPCS
        ----------------------------------------
        -0.0325040608617000  -0.0333250879034000
        -0.0470001557641000  -0.0458749339862000
        -0.0370348732066000  -0.0364138324117000
        -0.0359768558469000  -0.0358236703894000
        -0.0356064894317000  -0.0354277722246000
        -0.0355250220872000  -0.0353312047875000
        -0.0356105862451000  -0.0354251625379000

        The reference value is taken as the average of the values
        for Chorin and IPCS on the finest mesh.
        """
        if t < self.T:
            return 0.0

        return -0.0355
'''

if __name__ == "__main__":
    from demo_main import demo_main
    demo_main(PipeAneurysm)