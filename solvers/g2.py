__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-01-01"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian Valen-Sendstad, 2009.

# This implementation differs from the reference implementation g2ref.py
# in several ways (see comments in g2ref.py). It can be verified against
# g2ref.py by setting the tolerance to 1e-2:
#
#   if converged(r, iter, "Fixed-point", tolerance=1e-2): break
#
# The error should then be 0.0062193. With the default tolerance
# (1e-6), the error should be 46.503 when running ./ns g2ref g2.

from solverbase import *
from g2cppcode import *
g2ref = False

class Solver(SolverBase):
    "G2 (stabilized cG(1)cG(1)) by Hoffman and Johnson."

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):

        # Get problem parameters
        mesh = problem.mesh

        # Set time step
        dt, t, t_range = problem.timestep(problem)

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 1)
        Q = FunctionSpace(mesh, "CG", 1)
        DG = FunctionSpace(mesh, "DG", 0)
        DGv = VectorFunctionSpace(mesh, "DG", 0)

        # Get initial and boundary conditions
        u0, p0 = problem.initial_conditions(V, Q)
        bcu, bcp = problem.boundary_conditions(V, Q, t)

        # Remove boundary stress term is problem is periodic
        if is_periodic(bcp):
            beta = Constant(0)
        else:
            beta = Constant(1)

        # FIXME: This should be moved to problem.boundary_conditions()
        pbar = problem.pressure_bc(Q)

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
        n  = FacetNormal(mesh)
        f  = problem.f

	# Stabilization parameters
	C1  = 4.0
	C2  = 2.0
        d1 = Expression(cppcode=cppcode_d1)
        d2 = Expression(cppcode=cppcode_d2)

        # Velocity system
        U = 0.5*(u0 + u)
        P = p1
        Fv = inner(v, u - u0)*dx + k*inner(v, grad(U)*W)*dx \
           + k*inner(epsilon(v), sigma(U, P, nu))*dx \
           - beta*k*nu*inner(v, grad(U).T*n)*ds + k*inner(v, pbar*n)*ds \
           - k*inner(v, f)*dx \
           + d1*k*inner(grad(v)*W, grad(U)*W)*dx + d2*k*div(v)*div(U)*dx
        av = lhs(Fv)
        Lv = rhs(Fv)

        if g2ref == True:
            # Pressure system
            ap = d1*inner(grad(q), grad(p))*dx
            Lp = - q*div(u1)*dx
        else:
            # Pressure system
            ap = inner(grad(q), grad(p))*dx
            Lp = -(1/d1)* q*div(u1)*dx

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

            # Get boundary conditions
            bcu, bcp = problem.boundary_conditions(V, Q, t)

            # Update stabilization parameters
	    cputime = time()
            d1.update(u0, problem.nu, dt, C1)
            d2.update(u0, problem.nu, dt, C2)

            # Solve nonlinear system by fixed-point iteration
            for iter in range(maxiter):

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
		cputime = time()
		bv = assemble(Lv)
		[bc.apply(Av, bv) for bc in bcu]
		cputime = time()
		solve(Av, u1.vector(), bv, "gmres", "ilu")

                # Compute projection of u1 onto piecewise constants (mean-value on each cell)
                bw = assemble(Lw)
		solve(Aw, W.vector(), bw, "gmres", "ilu")

                # Reassemble velocity system
		cputime = time()
                Av = assemble(av)
                bv = assemble(Lv)
                [bc.apply(Av, bv) for bc in bcu]

                # Reassemble pressure system
                bp = assemble(Lp)
                if len(bcp) == 0: normalize(bp)
                [bc.apply(Ap, bp) for bc in bcp]

                # Compute residuals
		cputime = time()
                rv = residual(Av, u1.vector(), bv)
		cputime = time()
                rp = residual(Ap, p1.vector(), bp)
                r = sqrt(rv**2 + rp**2)

                # Check for convergence
                if g2ref == True:
                    if has_converged(r, iter, "Fixed-point", tolerance=1e-2): break
                else:
                    if has_converged(r, iter, "Fixed-point"): break
                if g2ref == True:
                    d1.update(u0, problem.nu, dt, C1)
                    d2.update(u0, problem.nu, dt, C2)

            # Update
            self.update(problem, t, u1, p1)
            if g2ref == True:
                print 'Norm of velocity vector:' , norm(u1.vector())
                print 'Norm of pressure vector:' , norm(p1.vector())
	    u0.assign(u1)

        return u1, p1

    def __str__(self):
        return "G2"
