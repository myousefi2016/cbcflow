__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2008-02-01"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2008-2010.

from solverbase import *

class Solver(SolverBase):
    "Incremental pressure-correction scheme."

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):

        if str(problem)=="Aneurysm":
            pc = "jacobi"
        else:
            pc = "ilu"

        # Get problem parameters
        mesh = problem.mesh
        dt, t, t_range = problem.timestep(problem)

        # Define function spaces
        V = FunctionSpace(mesh, "CG", 1)
        Q = FunctionSpace(mesh, "CG", 1)
        DG = FunctionSpace(mesh, "DG", 0)

        segregated = True

        # Get initial and boundary conditions
        ux0, uy0, uz0, p0 = problem.initial_conditions(V, Q, segregated)
        bcux, bcuy, bcuz, bcp = problem.boundary_conditions(V, Q, t, segregated)

        # Remove boundary stress term is problem is periodic
        if is_periodic(bcp):
            beta = Constant(0)
        else:
            beta = Constant(1)

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        ux = uy = uz = u 

        # Functions
        ux0 = interpolate(ux0, V)
        uy0 = interpolate(uy0, V)
        uz0 = interpolate(uz0, V)

        ux1 = Function(V)
        uy1 = Function(V)
        uz1 = Function(V)
        p0 = interpolate(p0, Q)
        p1 = interpolate(p0, Q)
        nu = Constant(problem.nu)
        k  = Constant(dt)
        f  = problem.f
        n  = FacetNormal(mesh)

        # Tentative velocity step
        Ux = 0.5*(ux0 + ux)
        Fx1 = (1/k)*inner(v, ux - ux0)*dx +  v*(ux0*ux0.dx(0) + uy0*ux0.dx(1) + uz0*ux0.dx(2))*dx \
#            + inner(epsilon(v), sigma(Ux, p0, nu))*dx \
            + inner(grad(v), grad(Ux))*dx \   # FIXME, include also pressure term
            + inner(v, p0*n[0])*ds \
            - v*f[0]*dx
        ax1 = lhs(Fx1)
        Lx1 = rhs(Fx1)

        Uy = 0.5*(uy0 + uy)
        Fy1 = (1/k)*inner(v, uy - uy0)*dx +  v*(ux0*uy0.dx(0) + uy0*uy0.dx(1) + uz0*uy0.dx(2))*dx \
#            + inner(epsilon(v), sigma(Uy, p0, nu))*dx \
            + inner(grad(v), grad(Uy))*dx \   # FIXME, include also pressure term
            + inner(v, p0*n[0])*ds \
            - v*f[1]*dx
        ay1 = lhs(Fy1)
        Ly1 = rhs(Fy1)

        Uz = 0.5*(uz0 + uz)
        Fz1 = (1/k)*inner(v, uz - uz0)*dx +  v*(ux0*uz0.dx(0) + uy0*uz0.dx(1) + uz0*uz0.dx(2))*dx \
#            + inner(epsilon(v), sigma(Uz, p0, nu))*dx \
            + inner(grad(v), grad(Uz))*dx \   # FIXME, include also pressure term
            + inner(v, p0*n[0])*ds \
            - v*f[2]*dx
        az1 = lhs(Fz1)
        Lz1 = rhs(Fz1)


        # Pressure correction
        a2 = inner(grad(q), grad(p))*dx
        L2 = inner(grad(q), grad(p0))*dx - (1/k)*q*(ux1.dx(0) + uy1.dx(1) + uz1.dx(2))*dx

        # Velocity correction
        a3 = inner(v, u)*dx
        L3 = inner(v, u1)*dx - k*inner(v, grad(p1 - p0))*dx

        # Assemble matrices
        # FIXME, in this case all matrices are equal!!!
        Ax1 = assemble(a1)
        Ay1 = assemble(a1)
        Az1 = assemble(a1)

        A2 = assemble(a2)
        A3 = assemble(a3)

        # Time loop
        self.start_timing()
        for t in t_range:

            # Get boundary conditions
            bcux, bcuy, bcuz, bcp = problem.boundary_conditions(V, Q, t)

            # Compute tentative velocity step
            bx = assemble(Lx1)
            by = assemble(Ly1)
            bz = assemble(Lz1)

            [bc.apply(Ax1, b) for bc in bcux]
            solve(Ax1, ux1.vector(), b, "gmres", "ilu")

            [bc.apply(Ay1, b) for bc in bcuy]
            solve(Ay1, uy1.vector(), b, "gmres", "ilu")

            [bc.apply(Az1, b) for bc in bcuz]
            solve(Az1, uz1.vector(), b, "gmres", "ilu")


            # Pressure correction
            b = assemble(L2)
            if len(bcp) == 0 or is_periodic(bcp): normalize(b)
            [bc.apply(A2, b) for bc in bcp]
            if is_periodic(bcp):
                solve(A2, p1.vector(), b)
            else:
                solve(A2, p1.vector(), b, 'gmres', 'hypre_amg')
            if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())

            # Velocity correction
            b = assemble(L3)
            [bc.apply(A3, b) for bc in bcu]
            solve(A3, u1.vector(), b, "gmres", pc)

            # Update
            self.update(problem, t, u1, p1)
            u0.assign(u1)
            p0.assign(p1)

        return u1, p1

    def __str__(self):
        return "IPCS"
