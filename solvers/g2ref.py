__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-01-01"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian Valen-Sendstad, 2009.

# This particular implementation of G2 is tailored to produce the
# exact same results as the Unicorn 3D benchmark problem found in
# unicorn-0.1.0/ucsolver/icns/bench2D3D for the cylinder.py test
# case.
#
# Run the test case by
#
#   ./ns g2ref g2ref
#
# The following reference values are obtained at T = 0.2:
#
# Norm of velocity vector: 106.5839601025451   (Unicorn, Kristian)
# Norm of velocity vector: 106.5839601165444   (Unicorn, Anders)
# Norm of velocity vector: 106.583960098       (g2ref.py, Kristian)
# Norm of velocity vector: 106.583960098       (g2ref.py, Anders)
#
# Norm of pressure vector: 38.25864714422254   (Unicorn, Kristian)
# Norm of pressure vector: 38.25864682055209   (Unicorn, Anders)
# Norm of pressure vector: 38.258646829        (g2ref.py, Kristian)
# Norm of pressure vector: 38.2586468266       (g2ref.py, Anders)
#
# New results (less strict default tolerance in new linear solvers in DOLFIN)
#
# Norm of velocity vector: 106.583962337
# Norm of pressure vector: 38.2586340271
#
# Differences from g2.py:
#
# 1. Different tolerance for iteration
# 2. Using FFC tensor representation
# 3. Optimized compilation
# 4. Not using div sigma formulation

# FIXME: Stabilization vector d1 is on LHS in pressure eq., but never reassembled. Should be on RHS which is reassembled.

from solverbase import *
from g2cppcode import *
from numpy import linspace

# Use same form compiler options as Unicorn to get comparable timings
parameters["form_compiler"]["representation"] = "tensor"
parameters["form_compiler"]["cpp_optimize"] = True

class Solver(SolverBase):
    "G2 (stabilized cG(1)cG(1)) by Hoffman and Johnson."

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):

        # Get problem parameters
        mesh = problem.mesh

        # FIXME: Is this correct Kristian?
        # Set time step
        dt = 0.021650635094
        t_range = linspace(0, 10*dt, 10)[1:]
        t = dt

        # Set parameters for method
        tol = 1.0e-2
        maxiter = 50

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 1)
        Q = FunctionSpace(mesh, "CG", 1)
        DG = FiniteElement("DG", "tetrahedron", 0)
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
        u0 = interpolate(u0, V)   # velocity at previous time step
        u1 = Function(V)          # current velocity
        p1 = interpolate(p0, Q)   # current pressure
        W  = Function(DGv)        # cell mean linearized velocity
        nu = Constant(problem.nu) # kinematic viscosity
        k  = Constant(dt)         # time step
        f  = problem.f            # body forces

	# Stabilization parameters
	C1  = 4.0
	C2  = 2.0
        d1 = Expression(cppcode_d1, element=DG)
        d2 = Expression(cppcode_d2, element=DG)

        # Update stabilization parameters
        d1.update(u0, problem.nu, dt, C1)
        d2.update(u0, problem.nu, dt, C2)

        # Velocity system
        U = 0.5*(u0 + u)
        Fv = inner(v, u - u0) + k*inner(v, grad(U)*W) + k*nu*inner(grad(v), grad(U))  - k*inner(div(v), p1) - k*inner(v, f) + \
             d1*k*inner(grad(v)*W, grad(U)*W) + d2*k*div(v)*div(U)
        av = lhs(Fv*dx)
        Lv = rhs(Fv*dx)

        # Pressure system
        ap = d1*inner(grad(q), grad(p))*dx
        Lp = - q*div(u1)*dx

        # Projection of velocity
        aw = inner(z, w)*dx
        Lw = inner(z, u1)*dx

        # Assemble matrices
        Av = assemble(av)
        Ap = assemble(ap)
        Aw = assemble(aw)

        # Time loop
        self.start_timing()
        for (time_step, t) in enumerate(t_range):
            print "============================================================================="
    	    print "Starting time step %d, t = %g and dt = %g" % (time_step, t, dt)
            print

            # Solve nonlinear system by fixed-point iteration
            for iter in range(maxiter):
                # Update stabilization parameters
                cputime = time()
                d1.update(u0, problem.nu, dt, C1)
                d2.update(u0, problem.nu, dt, C2)
                print "Computed stabilization in:", time() - cputime
                print

                print "-----------------------------------------------------------------------------"
                print "Starting iteration", iter
                print

                # Compute pressure
		cputime = time()
                bp = assemble(Lp)
                if len(bcp) == 0: normalize(bp)
                [bc.apply(Ap, bp) for bc in bcp]
		print "Assembled pressure RHS in:", time() - cputime
		cputime = time()
                solve(Ap, p1.vector(), bp, "gmres", "amg_hypre")
                if len(bcp) == 0: normalize(p1.vector())
		print "Solved pressure in: ", time() - cputime
		print "Pressure norm =", norm(p1.vector())
		print

                # Compute velocity
		cputime = time()
		bv = assemble(Lv)
		[bc.apply(Av, bv) for bc in bcu]
		print "Assembled velocity RHS in:", time() - cputime
		cputime = time()
		solve(Av, u1.vector(), bv, "gmres", "ilu")
		print "Solved velocity in:", time() - cputime
		print "Velocity norm =", norm(u1.vector())
                print

                # Compute projection of u1 onto piecewise constants (mean-value on each cell)
                bw = assemble(Lw)
		solve(Aw, W.vector(), bw, "gmres", "ilu")

                # Reassemble velocity system
		cputime = time()
                Av = assemble(av)
                bv = assemble(Lv)
                [bc.apply(Av, bv) for bc in bcu]
		print "Assembled velocity LSH and RHS in:" , time() - cputime

                # Reassemble pressure system
                bp = assemble(Lp)
                if len(bcp) == 0: normalize(bp)
                [bc.apply(Ap, bp) for bc in bcp]

                # Compute residuals
		cputime = time()
                rv = residual(Av, u1.vector(), bv)
		print "Computed residual for velocity in:" , time() - cputime
		cputime = time()
                rp = residual(Ap, p1.vector(), bp)
		print "Computed residual for pressure in:" , time() - cputime
                r = sqrt(rv*rv + rp*rp)
                print

                # Check convergence
                print "Residual =", r
		print ""

                if r < tol: break
                if iter == maxiter - 1: raise RuntimeError, "Fixed-point iteration did not converge."

	    print 'Norm of velocity vector:' , norm(u1.vector())
	    print 'Norm of pressure vector:' , norm(p1.vector())
            print

            # Update
            self.update(problem, t, u1, p1)
	    u0.assign(u1)

        return u1, p1

    def __str__(self):
        return "G2ref"
