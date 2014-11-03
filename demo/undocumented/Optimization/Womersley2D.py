#!/usr/bin/env python

from cbcflow import *
from cbcflow.dol import *
from cbcpost import ParamDict, PostProcessor

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
    "2D pipe test problem with known stationary analytical solution."

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
        self.alpha = 2.0 * Q / (pi * RADIUS**4)
        self.beta = 2.0 * self.nu * self.alpha

        # Toggle between constant and transient flow rate
        if 0:
            self.Q_coeffs = [(0.0, Q), (1.0, Q)]
        else:
            T = self.params.T
            P = self.params.period
            tvalues = np.linspace(0.0, P)
            Qfloor = 0.25
            Qpeak = 1.0
            Qvalues = Q * (Qfloor + (Qpeak-Qfloor)*np.sin(pi*np.mod((P-tvalues)/P, P)**3)**2)
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
            rho = 1.0,
            mu=1.0/30.0,
            )
        params.update(
            # Spatial parameters
            refinement_level=2,
            # Analytical solution parameters
            Q=1.0,
            coeffstype="Q",
            num_womersley_coefficients=25,
            )
        return params

    def womersley_solution(self, spaces, t):
        # Create womersley objects
        ua = make_womersley_bcs(self.Q_coeffs, self.mesh, self.left_boundary_id, self.nu, None, self.facet_domains,
                                self.params.coeffstype, num_fourier_coefficients=self.params.num_womersley_coefficients)
        for uc in ua:
            uc.set_t(t)

        pa = Expression("-beta * x[0]", beta=1.0)
        pa.beta = self.beta # TODO: This is not correct unless stationary...
        return (ua, pa)

    def __poiseuille_solution(self, spaces, t):
        A = 2*RADIUS
        Q = self.params.Q
        mu = self.params.mu

        dpdx = 3*Q*mu/(A*RADIUS**2)

        ux = Expression("(radius*radius - x[1]*x[1])/(2*mu)*dpdx", radius=RADIUS, mu=mu, dpdx=dpdx, degree=2)
        uy = Constant(0.0)

        u = [ux, uy]

        p = Expression("dpdx * (length-x[0])", dpdx=dpdx, length=LENGTH)

        return (u, p)

    def observations(self, spaces, t):
        return []

    def controls(self, spaces):
        return []

    def cost_functionals(self, spaces, t, observations, controls):
        return []

    def initial_conditions(self, spaces, controls):
        d = spaces.d

        # TODO: Get icu, icp from controls
        icu = [Function(spaces.U) for i in range(d)]
        icp = Function(spaces.Q)

        # Interpolate womersley solution into entire initial condition (only possible for this simple pipe case!)
        ua, pa = self.womersley_solution(spaces, 0.0)
        for i in range(d):
            icu[i].interpolate(ua[i])

        return icu, icp

    def boundary_conditions(self, spaces, u, p, t, controls):
        # Create no-slip bcs
        d = len(u)
        u0 = [Constant(0.0) for i in range(d)]
        noslip = (u0, self.wall_boundary_id, "strong")

        # Create inflow bcs
        uin = [Function(spaces.U) for i in range(d)] # TODO: Get uin from controls
        inflow = (uin, self.left_boundary_id, "nietche")

        # Interpolate Womersley profile into inflow bc functions
        self.ua, self.pa = self.womersley_solution(spaces, t)
        for iucomp, ucomp in zip(uin, self.ua):
            iucomp.interpolate(ucomp)

        # Create outflow bcs for pressure
        outflow = (Constant(0.0), self.right_boundary_id, "natural")

        # Return bcs in two lists
        bcu = [inflow, noslip]
        bcp = [outflow]
        return (bcu, bcp)

    def update(self, spaces, u, p, t, timestep, bcs, observations, controls, cost_functionals):
        # TODO: Drop this for optimization run
        # Update time in boundary condition expressions
        bcu, bcp = bcs
        inflow, noslip = bcu
        uin = inflow[0]
        for ucomp in self.ua:
            ucomp.set_t(t)
        for iucomp, ucomp in zip(uin, self.ua):
            iucomp.interpolate(ucomp)

        # TODO: Update observations

        # TODO: Update initial guess for controls at time t
        if controls:
            m0 = controls[0]
            V = m0.function_space()
            m = Function(V)
            #m.interpolate(initial_control_guess)
            controls.append(m)
            m0.assign(m)

        # TODO: Add contribution to cost functionals at time t

def main():
    set_log_level(100)

    problem = Womersley2D(
        ParamDict(
            dt=1e-3,
            T=0.2,#8,
            num_periods=None,
            refinement_level=2,
            )
        )

    scheme = CoupledPicard(
        ParamDict(
            nietche=ParamDict(
                enable=True,
                formulation=1,
                stabilize=True,
                gamma=100.0,
                ),
            scale_by_dt=True,
            enable_convection=True, # False = Stokes
            )
        )
    params_string = "__".join("{}_{}".format(k, scheme.params.nietche[k]) for k in scheme.params.nietche)
    equation = "NavierStokes" if scheme.params.enable_convection else "Stokes"
    casedir = "results_demo_%s_%s_%s_%s" % (problem.shortname(), scheme.shortname(), equation, params_string)
    plot_and_save = dict(plot=True, save=True)
    fields = [
        Pressure(plot_and_save),
        Velocity(plot_and_save),
        LocalCfl(plot_and_save),
        ]
    postproc = PostProcessor({"casedir": casedir})
    postproc.add_fields(fields)

    solver = NSSolver(problem, scheme, postproc)
    solver.solve()


if __name__ == "__main__":
    main()
