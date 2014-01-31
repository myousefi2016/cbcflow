#!/usr/bin/env python
__author__ = "Martin Sandve Alnes <martinal@simula.no>"
__date__ = "2013-08-12"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from cbcflow import *
from cbcflow.dol import *

import numpy as np

LENGTH = 10.0
RADIUS = 0.5

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > LENGTH*(1.0 - DOLFIN_EPS)

class Womersley2D(NSProblem):
    "2D pipe test problem with known transient analytical solution."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Create mesh
        refinements = [4,8,16,32,64]
        N = refinements[self.params.refinement_level]
        M = int(N*LENGTH/(2*RADIUS) + 0.5)
        mesh = UnitSquareMesh(M, N)
        x = mesh.coordinates()[:,0]
        y = mesh.coordinates()[:,1]
        x = LENGTH*x
        y = RADIUS*2*(y - 0.5)
        mesh.coordinates()[:,0] = x
        mesh.coordinates()[:,1] = y
        
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
        Q = self.params.Q
        self.nu = self.params.mu / self.params.rho

        # Beta is the Poiseuille pressure drop if the flow rate is stationary Q
        self.beta = 4.0 * self.nu * Q / (pi * RADIUS**4)

        # Setup transient flow rate coefficients
        if 0:
            print "Using stationary bcs."
            self.Q_coeffs = [(0.0, Q), (1.0, Q)]
        else:
            print "Using transient bcs."
            T = self.params.T
            P = self.params.period
            tvalues = np.linspace(0.0, P)
            #Qfloor, Qpeak = 1.0, 0.0
            #Qfloor, Qpeak = 0.3, 0.7
            Qfloor, Qpeak = -0.2, 1.0
            #Qfloor, Qpeak = 1.0, 1.0
            Qvalues = Q * (Qfloor + (Qpeak-Qfloor)*np.sin(pi*((P-tvalues)/P)**2)**2)
            self.Q_coeffs = zip(tvalues, Qvalues)

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

        # To be able to compare womersley profiles based on midpoint velocity and flow rate,
        # we here sample the flow rate based womersley profile in the midpoints to get
        # matching midpoint velocity coefficients
        if self.params.coeffstype == "V":
            # Build midpoint velocity coefficients from evaluating flow rate based womersley profile
            ua = make_womersley_bcs(self.Q_coeffs, self.mesh, self.left_boundary_id, self.nu, None, self.facet_domains, "Q",
                                    num_fourier_coefficients=self.params.num_womersley_coefficients)
            ua = ua[0]
            V_coeffs = []
            x = ua.center
            value = np.zeros((1,))
            for t,Q in self.Q_coeffs:
                ua.set_t(t)
                ua.eval(value, x)
                V_coeffs.append((float(t), float(value[0])))
            self.Q_coeffs = V_coeffs

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=None,
            dt=1e-2,
            period=0.8,
            num_periods=1.0,
            # Physical parameters
            rho=1.0,
            mu=1.0/30.0,
           #mu=1.0e6, #1.0/30.0,
            )
        params.update(
            # Spatial parameters
            #N=32,
            refinement_level=0,
            # Analytical solution parameters
            Q=1.0,
            coeffstype="Q",
            num_womersley_coefficients=25,
            )
        return params

    def analytical_solution(self, spaces, t):
        # Create womersley objects
        ua = make_womersley_bcs(self.Q_coeffs, self.mesh, self.left_boundary_id, self.nu, None, self.facet_domains,
                                self.params.coeffstype, num_fourier_coefficients=self.params.num_womersley_coefficients)
        #ua = womersley(self.Q_coeffs, self.mesh, self.facet_domains, self.left_boundary_id, self.nu) # TODO
        for uc in ua:
            uc.set_t(t)           
            
        pa = Expression("-beta * x[0]", beta=1.0)
        pa.beta = self.beta # TODO: This is not correct unless stationary...
        return (ua, pa)
    
    def test_fields(self):
        return [Velocity(), Pressure()]
    
    def test_references(self, spaces, t):
        # Create womersley objects
        ua = make_womersley_bcs(self.Q_coeffs, self.mesh, self.left_boundary_id, self.nu, None, self.facet_domains,
                                self.params.coeffstype, num_fourier_coefficients=self.params.num_womersley_coefficients)
        #ua = womersley(self.Q_coeffs, self.mesh, self.facet_domains, self.left_boundary_id, self.nu) # TODO
        for uc in ua:
            uc.set_t(t)
        pa = Expression("-beta * x[0]", beta=1.0)
        pa.beta = self.beta # TODO: This is not correct unless stationary...

        return [as_vector(ua), pa]

    def initial_conditions(self, spaces, controls):
        return self.analytical_solution(spaces, 0.0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        # Create no-slip bcs
        d = len(u)
        u0 = [Constant(0.0)] * d
        noslip = (u0, self.wall_boundary_id)

        # Get other bcs from analytical solution functions
        ua, pa = self.analytical_solution(spaces, t)

        # Create inflow boundary conditions for velocity
        inflow = (ua, self.left_boundary_id)

        # Create outflow boundary conditions for velocity
        u_outflow = (ua, self.right_boundary_id)

        # Create outflow boundary conditions for pressure
        p_outflow = (pa, self.right_boundary_id)

        # Return bcs in two lists
        bcu = [noslip, inflow]
        bcp = []
        if 1: # Switch between pressure or dbc at outlet
            bcp += [p_outflow]
        else:
            bcu += [u_outflow]
        return (bcu, bcp)

    def update(self, spaces, u, p, t, timestep, bcs, observations, controls):
        bcu, bcp = bcs
        uin = bcu[1][0]
        for ucomp in uin:
            ucomp.set_t(t)


def main():
    problem = Womersley2D()
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
