__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-01-01"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian Valen-Sendstad, 2009.

# This file contains two different versions of the G2 method.
#
# 1. One version which is tailored to produce the exact same results
# as Unicorn 0.1.0 for the 3D benchmark problem available in the
# directory unicorn-0.1.0/ucsolver/icns/bench2D3D.
#
# 2. Another version which differs from the reference implementation
# in a number of ways: stricter tolerance, div-sigma formulation,
# divide by delta_1 in pressure system (so it gets included in RHS and
# is thereby reassembled), and updating the stabilization parameters
# outside of the iteration loop (for efficiency).
#
# One may select either of the two methods by setting the flag 'g2ref'
# to either True or False (default). The flag is automatically set to
# True for the 'g2ref' test problem.
#
# Reference results obtained with Unicorn 0.1.0:
#
# Norm of velocity vector: 106.5839601165444  (Unicorn, Anders)
# Norm of velocity vector: 106.5839601025451  (Unicorn, Kristian)
# Norm of pressure vector: 38.25864682055209  (Unicorn, Anders)
# Norm of pressure vector: 38.25864714422254  (Unicorn, Kristian)

from solverbase import *
from g2cppcode import *

class Solver(SolverBase):
    "G2 (stabilized cG(1)cG(1)) by Hoffman and Johnson."

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):

        # Check if we should use reference implementation (Unicorn)
        g2ref = str(problem) == "G2ref"

        # Get problem parameters
        mesh = problem.mesh

        # Set time step
        dt, t, t_range = problem.timestep(problem)

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 1)
        Q = FunctionSpace(mesh, "CG", 1)
        DG = FiniteElement("DG", None, 0)
        DGv = VectorFunctionSpace(mesh, "DG", 0)

        # Get initial and boundary conditions
        u0, p0 = problem.initial_conditions(V, Q)
        bcu, bcp = problem.boundary_conditions(V, Q, t)

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        z = TestFunction(DGv)
        u = TrialFunction(V)
        p = TrialFunction(Q)
        w = TrialFunction(DGv)

        # Functions
        u0 = interpolate(u0, V)
        u1 = Function(V)
        p1 = interpolate(p0, Q)
        W  = Function(DGv)
        nu = Constant(problem.nu)
        k  = Constant(dt)
        f  = problem.f
        n  = FacetNormal(mesh)

	# Stabilization parameters
	C1  = 4.0
	C2  = 2.0
        d1 = Expression(cppcode_d1, element=DG)
        d2 = Expression(cppcode_d2, element=DG)

        # Remove boundary stress term if problem is periodic
        if is_periodic(bcp):
            beta = Constant(0)
        else:
            beta = Constant(1)

        # FIXME: This should be moved to problem.boundary_conditions()
        pbar = problem.pressure_bc(Q)

        # Velocity system
        U = 0.5*(u0 + u)
        P = p1
        if g2ref:
            stress_terms = k*nu*inner(grad(v), grad(U))*dx - k*inner(div(v), P)*dx
        else:
            stress_terms = k*inner(epsilon(v), sigma(U, P, nu))*dx \
                         - beta*k*nu*inner(v, grad(U).T*n)*ds + k*inner(v, pbar*n)*ds
        Fv = inner(v, u - u0)*dx + k*inner(v, grad(U)*W)*dx \
           + stress_terms - k*inner(v, f)*dx \
           + d1*k*inner(grad(v)*W, grad(U)*W)*dx + d2*k*div(v)*div(U)*dx
        av = lhs(Fv)
        Lv = rhs(Fv)

        # Pressure system
        if g2ref:
            ap = d1*inner(grad(q), grad(p))*dx
            Lp = -q*div(u1)*dx
        else:
            ap = inner(grad(q), grad(p))*dx
            Lp = -(1/d1)*q*div(u1)*dx

        # Projection of velocity
        aw = inner(z, w)*dx
        Lw = inner(z, u1)*dx

        # Update stabilization parameters
        d1.update(u0, problem.nu, dt, C1)
        d2.update(u0, problem.nu, dt, C2)

        # Assemble matrices
        Av = assemble(av)
        Ap = assemble(ap)
        Aw = assemble(aw)

        # Time loop
        self.start_timing()
        for t in t_range:

            # FIXME: Update boundary conditions here in all solvers
            # Update boundary conditions
            bcu, bcp = problem.boundary_conditions(V, Q, t)

            # Update stabilization parameters
            if not g2ref:
                d1.update(u0, problem.nu, dt, C1)
                d2.update(u0, problem.nu, dt, C2)

            # Solve nonlinear system by fixed-point iteration
            for iter in range(maxiter):

                # Update stabilization parameters
                if g2ref:
                    d1.update(u0, problem.nu, dt, C1)
                    d2.update(u0, problem.nu, dt, C2)

                # Compute pressure
                bp = assemble(Lp)
                if len(bcp) == 0 or is_periodic(bcp): normalize(bp)
                [bc.apply(Ap, bp) for bc in bcp]
                if is_periodic(bcp):
                    solve(Ap, p1.vector(), bp)
                else:
                    solve(Ap, p1.vector(), bp, "gmres", "amg_hypre")
                if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())

                # Compute velocity
		bv = assemble(Lv)
		[bc.apply(Av, bv) for bc in bcu]
		solve(Av, u1.vector(), bv, "gmres", "ilu")

                # Compute projection of u1 onto piecewise constants (mean-value on each cell)
                bw = assemble(Lw)
		solve(Aw, W.vector(), bw, "gmres", "ilu")

                # Reassemble velocity system
                Av = assemble(av)
                bv = assemble(Lv)
                [bc.apply(Av, bv) for bc in bcu]

                # Reassemble pressure system
                bp = assemble(Lp)
                if len(bcp) == 0: normalize(bp)
                [bc.apply(Ap, bp) for bc in bcp]

                # Compute residuals
                rv = residual(Av, u1.vector(), bv)
                rp = residual(Ap, p1.vector(), bp)
                r = sqrt(rv**2 + rp**2)

                # Check for convergence
                if g2ref:
                    if has_converged(r, iter, "Fixed-point", tolerance=1e-2): break
                else:
                    if has_converged(r, iter, "Fixed-point"): break

            # Update
            self.update(problem, t, u1, p1)
	    u0.assign(u1)

        # Compare solution with reference implementation (Unicorn)
        if g2ref:
            check_g2(u1, p1)

        return u1, p1

    def __str__(self):
        return "G2ref"

def check_g2(u, p):
    "Compare solution with reference implementation (Unicorn)."
    unorm = norm(u.vector())
    pnorm = norm(p.vector())
    unorm_ref = 106.5839601165444
    pnorm_ref = 38.25864682055209
    error = sqrt((unorm - unorm_ref)**2 + (pnorm - pnorm_ref)**2)
    print "Norm of velocity vector: %.16g (should be %.16g)" % (unorm, unorm_ref)
    print "Norm of pressure vector: %.16g (should be %.16g)" % (pnorm, pnorm_ref)
    print "Difference from reference implementation: %g" % error
    if error < 1e-7:
        print "G2 looks OK!"
    else:
        print "*** ERROR: Something is wrong with G2"
