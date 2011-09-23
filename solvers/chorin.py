__author__ = "Anders Logg <logg@simula.no> and Kent-Andre Mardal <kent-and@simula.no>"
__date__ = "2008-03-03"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class Solver(SolverBase):
    "Original pressure-correction scheme by Chorin and Temam."

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):

        # Get problem parameters
        mesh = problem.mesh
        dt, t, t_range = problem.timestep(problem)

        if str(problem)=="Aneurysm":
            pc = "jacobi"
        else:
            pc = "ilu"

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)
        DG = FunctionSpace(mesh, "DG", 0)

        # Get initial and boundary conditions
        u0, p0 = problem.initial_conditions(V, Q)
        bcu, bcp = problem.boundary_conditions(V, Q, t)

        # Define test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # Define functions
        us = Function(V)
        u0 = interpolate(u0, V)
        u1 = interpolate(u0, V)
        p1 = interpolate(p0, Q)
        nu = Constant(problem.nu)
        k  = Constant(dt)
        f  = problem.f

        # Tentative velocity step
        F1 = (1/k)*inner(v, u - u0)*dx + inner(v, grad(u0)*u0)*dx \
            + nu*inner(grad(v), grad(u))*dx - inner(v, f)*dx
        a1 = lhs(F1)
        L1 = rhs(F1)

        # Poisson problem for the pressure
        a2 = inner(grad(q), grad(p))*dx
        L2 = -(1/k)*q*div(us)*dx

        # Velocity update
        a3 = inner(v, u)*dx
        L3 = inner(v, us)*dx - k*inner(v, grad(p1))*dx

        # Assemble matrices
        A1 = assemble(a1)
        A2 = assemble(a2)
        A3 = assemble(a3)

        # Time loop
        self.start_timing()
        for t in t_range:

            # Compute tentative velocity
            bcu, bcp = problem.boundary_conditions(V, Q, t)
            b = assemble(L1)
            [bc.apply(A1, b) for bc in bcu]
            solve(A1, us.vector(), b, 'gmres', 'ilu')

            # Compute p1
            b = assemble(L2)
            if len(bcp) == 0 or is_periodic(bcp): normalize(b)
            [bc.apply(A2, b) for bc in bcp]
            if is_periodic(bcp):
                solve(A2, p1.vector(), b)
            else:
                solve(A2, p1.vector(), b, 'cg', 'hypre_amg')
            if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())

            # Compute u1
            b = assemble(L3)
            [bc.apply(A3, b) for bc in bcu]
            solve(A3, u1.vector(), b, 'gmres', pc)

            # Update
            self.update(problem, t, u1, p1)
            u0.assign(u1)

        return u1, p1

    def __str__(self):
        return "Chorin"
