# Copyright (C) 2010-2014 Simula Research Laboratory
#
# This file is part of CBCFLOW.
#
# CBCFLOW is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CBCFLOW is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CBCFLOW. If not, see <http://www.gnu.org/licenses/>.
r"""
This scheme follows the same logic as in :class:`.IPCS`, but with a few notable exceptions.

A parameter :math:`\theta` is added to the diffusion and convection terms,
allowing for different evaluation of these, and the convection is handled semi-implicitly:

.. math::
    \frac{1}{\Delta t}\left( \tilde{u}^{n+1}-u^{n} \right)-
    \nabla\cdot\nu\nabla \tilde{u}^{n+\theta}+
    u^*\cdot\nabla \tilde{u}^{n+\theta}+\nabla p^{n}=f^{n+1},

where

.. math::
    u^* = \frac{3}{2}u^n - \frac{1}{2}u^{n-1}, \\
    \tilde{u}^{n+\theta} = \theta \tilde{u}^{n+1}+\left(1-\theta\right)u^n.

This convection term is unconditionally stable, and with :math:`\theta=0.5`,
this equation is second order in time and space [1]_.


In addition, the solution process is significantly faster by solving for each of the
velocity components separately, making for D number of smaller linear systems compared
to a large system D times the size.



.. [1] Simo, J. C., and F. Armero. *Unconditional stability and long-term behavior
    of transient algorithms for the incompressible Navier-Stokes and Euler equations.*
    Computer Methods in Applied Mechanics and Engineering 111.1 (1994): 111-154.

"""

from __future__ import division

from cbcpost.utils import cbc_log
from cbcpost.utils import get_memory_usage
from cbcflow.core.nsscheme import *

from cbcflow.schemes.utils import (compute_regular_timesteps,
                                   assign_ics_segregated,
                                   make_segregated_velocity_bcs,
                                   make_pressure_bcs,
                                   NSSpacePoolSegregated,
                                   RhsGenerator)
import petsc4py
import sys
from fenicstools.WeightedGradient import compiled_gradient_module, weighted_gradient_matrix
def delta_memory():
    if not hasattr(delta_memory, "M"):
        delta_memory.M = MPI.sum(mpi_comm_world(), get_memory_usage())
    
    Mtot = MPI.sum(mpi_comm_world(), get_memory_usage())
    delta = Mtot - delta_memory.M
    delta_memory.M = Mtot
    
    #return delta, Mtot
    return Mtot, delta

def _get_weighted_gradient(mesh, dims, v,p):
    DG = FunctionSpace(mesh, "DG", 0)
    q = TestFunction(DG)
    #print "WG: created DG", delta_memory()
    #G = assemble(TrialFunction(DG)*v*dx)
    #A = G
    A = assemble(TrialFunction(DG)*v*dx)
    dg = Function(DG)
    
    dPdX = []
    for d in dims:
        dP = assemble(p.dx(d)*q*dx)
        dPmat = as_backend_type(dP).mat()
        compiled_gradient_module.compute_DG0_to_CG_weight_matrix(A, dg)
        Amat = as_backend_type(A).mat()    

        Cpmat = Amat.matMultSymbolic(dPmat)
        Amat.matMultNumeric(dPmat, Cpmat)
        #print "WG: matMult %d" %d, delta_memory()

        # Perform some strange copies that apparently saves memory
        Cpmat2 = Cpmat.copy()
        Cpmat.destroy()
        Cp = PETScMatrix(Cpmat2)
        Cp = PETScMatrix(Cp)

        dPdX.append(Cp)

        MPI.barrier(mpi_comm_world())
        # Destroy petsc4py objects
        dPmat.destroy()
        Amat.destroy()
        Cpmat2.destroy()

    #del DG
    #del dg
    #del G
    #del A  

    return dPdX

class IPCS_Stable2(NSScheme):
    "Incremental pressure-correction scheme, fast and stable version."

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.update(
            # Default to P1-P1
            u_degree = 1,
            p_degree = 1,
            theta = 0.5,
            rebuild_prec_frequency = 1e16,
            u_tent_prec_structure = "same_nonzero_pattern",
            u_tent_solver_parameters = {},
            p_corr_solver_parameters = {},
            u_corr_solver_parameters = {},
            low_memory_version = False,
            store_rhs_matrix_p_corr = True,
            store_rhs_matrix_u_tent = True,
            store_rhs_matrix_u_corr = True,
            store_stiffness_matrix = True,
            #assemble_convection = "standard", # unassembled, debug
            )
        return params

    def solve(self, problem, timer):
        
        print "Begin solve: ", delta_memory()
        # Get problem parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n  = FacetNormal(mesh)
        dims = range(mesh.topology().dim())
        theta = self.params.theta

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Constant(timesteps[start_timestep], name="TIME")

        # Define function spaces
        spaces = NSSpacePoolSegregated(mesh, self.params.u_degree, self.params.p_degree)
        U = spaces.U
        Q = spaces.Q
        print "Created function spaces: ", delta_memory()

        # Test and trial functions
        v = TestFunction(U)
        q = TestFunction(Q)
        u = TrialFunction(U)
        p = TrialFunction(Q)

        # Functions
        u0 = as_vector([Function(U, name="u0_%d"%d) for d in dims]) # u^n
        u1 = as_vector([Function(U, name="u1_%d"%d) for d in dims]) # u^{n+1}
        u_ab = as_vector([Function(U, name="u_ab_%d"%d) for d in dims]) # Adams-Bashforth convection

        p0 = Function(Q, name="p0")
        p1 = Function(Q, name="p1")
        print "Created functions: ", delta_memory()

        # Get functions for data assimilation
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)

        # Apply initial conditions and use it as initial guess
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_segregated(u0, p0, spaces, ics)
        for d in dims:
            u1[d].assign(u0[d])

        # Update Adams-Bashford term for first timestep
        for d in dims:
            u_ab[d].vector().zero()
            u_ab[d].vector().axpy(1.5, u1[d].vector())
            u_ab[d].vector().axpy(-0.5, u0[d].vector())

        #for d in dims: u2[d].assign(u1[d])
        p1.assign(p0)

        # Make scheme-specific representation of bcs
        bcs = problem.boundary_conditions(spaces, u1, p1, t, controls)
        bcu = make_segregated_velocity_bcs(problem, spaces, bcs)
        bcp = make_pressure_bcs(problem, spaces, bcs)

        # Problem coefficients
        nu = Constant(problem.params.mu/problem.params.rho)
        rho = float(problem.params.rho)
        k  = Constant(dt)
        f  = as_vector(problem.body_force(spaces, t))

        print "done creating function stuff: ", delta_memory()
        timer.completed("create function spaces, functions and boundary conditions")

        # Tentative velocity step. Crank-Nicholson time-stepping is used for diffusion and convection.
        #a1 = (1/k) * inner(v, u) * dx()
        #a2 = inner(grad(v), nu*grad(u)) * dx()

        # Convection linearized as in Simo/Armero (1994)
        # Will set u_ab = 1.5*u1[r] - 0.5*u0[r]
        #a_conv = v * sum(u_ab[r] * u.dx(r) for r in dims) * dx()
        a1 = inner(u,v)*dx()
        a2 = inner(grad(v), nu*grad(u))*dx()
        a_conv  = inner(v, dot(u_ab, nabla_grad(u)))*dx()
        """
        if theta < 1.0:
            # Set a_conv to match rhs theta-weighting for RHSGenerator
            a_conv = Constant(1-theta)*a_conv
            Kconv_axpy_factor = theta/(1-theta)
        else:
            Kconv_axpy_factor = 1.0

        Kconv = Matrix() # assembled from a_conv in the time loop
        """
        
        #M = assemble(inner(u,v)*dx())


        print "before A_u_tent assemble", delta_memory()
        A_u_tent = assemble(1/k*inner(u,v)*dx())
        print "after A_u_tent assemble", delta_memory()
        
        #print "before A_conv+A assemble", MPI.sum(mpi_comm_world(), get_memory_usage())
        if not self.params.low_memory_version and self.params.store_stiffness_matrix:
            A = assemble(a2)
            K_conv = assemble(a_conv)
        else:
            print "Will not store stiffness matrix"
            A = assemble(a2+a_conv)

        print "after A_conv+A assemble", delta_memory()
        #A_u_tent = Matrix(M+A+K_conv)
        #A_u_tent = M+A+K_conv
        #A_u_tent.zero()
        
        #A_u_tent = assemble(M)
        #import ipdb; ipdb.set_trace()
        #A_u_tent.axpy(1.0, M, False)
        #A_u_tent.axpy(1.0, A, False)
        #A_u_tent.axpy(1.0, K_conv, False)
        
        #for bc in bcu:
        #    M *= 1/dt
        #    bc[0].apply(A_u_tent)
        #    M *= dt
        #    bc[0].apply(A)
        #A_u_tent.axpy(1.0, M, False)
        # Create the static part of the coefficient matrix for the tentative
        # velocity step. The convection is added in the time loop. We use a
        # single matrix for all dimensions, which means that the BCs must match
        # (they must apply to the full vector space at the vertex, not a
        # subspace.
        """
        print "before A_u_tent matrix creation: ", MPI.sum(mpi_comm_world(), get_memory_usage())
        A_u_tent = assemble(a1+theta*a2)
        print "after A_u_tent matrix creation: ", MPI.sum(mpi_comm_world(), get_memory_usage())
        print "after first matrix creation: ", MPI.sum(mpi_comm_world(), get_memory_usage())
        # Create matrices for generating the RHS
        print "before B matrix creation: ", MPI.sum(mpi_comm_world(), get_memory_usage())
        B = assemble(a1-(1-theta)*a2)
        print "after B matrix creation: ", MPI.sum(mpi_comm_world(), get_memory_usage())
        print "before M matrix creation: ", MPI.sum(mpi_comm_world(), get_memory_usage())
        M = assemble(v*u*dx())
        print "before M matrix creation: ", MPI.sum(mpi_comm_world(), get_memory_usage())
        """
        # Define how to create the RHS for the tentative velocity. The RHS is
        # (of course) different for each dimension.
        
        rhs_u_tent = [None]*len(dims)
        if not self.params.low_memory_version and self.params.store_rhs_matrix_u_tent:
        #if True:
            for d in dims:
                print "before C creation: ", delta_memory()
                C = assemble(-v*p*n[d]*ds() + v.dx(d)*p*dx())
                print "after C creation: ", delta_memory()
                #C = assemble(-v*p.dx(d)*dx())
                rhs_u_tent[d] = RhsGenerator(U)
                #rhs_u_tent[d] += B, u0[d]
                rhs_u_tent[d] += A_u_tent, u0[d]
                rhs_u_tent[d] += C, p0
                #rhs_u_tent[d] += M, f[d]
                #if theta < 1.0:
                #    rhs_u_tent[d] -= Kconv, u0[d]
        else:
            print "Will not store rhs matrix for u_tent"
            #b_tent = [Vector() for d in dims]
            rhs_u_tent = [lambda: A_u_tent*u0[d].vector()+assemble(-v*p0*n[d]*ds() + v.dx(d)*p0*dx()) for d in dims]

        # Apply BCs to LHS
        #for bc in bcu:
        #    bc[0].apply(A_u_tent)

        # Tentative velocity solver
        solver_u_tent = LinearSolver(*self.params.solver_u_tent)
        #solver_u_tent = PETScKrylovSolver(*self.params.solver_u_tent)
        #solver_u_tent.parameters["monitor_convergence"] = True
        #solver_u_tent.set_operator(A_u_tent)
        if 'preconditioner' in solver_u_tent.parameters:
                solver_u_tent.parameters['preconditioner']['structure'] = 'same'
        solver_u_tent.parameters.update(self.params.u_tent_solver_parameters)
        print "create tenative velocity solver", delta_memory()
        timer.completed("create tenative velocity solver")

        # Pressure correction
        A_p_corr = assemble(inner(grad(q), grad(p))*dx())
        if not self.params.low_memory_version and self.params.store_rhs_matrix_p_corr:
            rhs_p_corr = RhsGenerator(Q)
            rhs_p_corr += A_p_corr, p0
            Ku = [None]*len(dims)
            for d in dims:
                print "before Ku creation: ", delta_memory()
                Ku[d] = assemble(-(1/k)*q*u.dx(d)*dx()) # TODO: Store forms in list, this is copied below
                print "after Ku creation: ", delta_memory()
                #print Ku[d].norm('frobenius')
                rhs_p_corr += Ku[d], u1[d]
        else:
            rhs_p_corr = lambda: A_p_corr*p0.vector() + assemble(-(1/k)*q*sum(u1[d].dx(d) for d in dims)*dx())

        # Pressure correction solver
        if self.params.solver_p:
            solver_p_params = self.params.solver_p
        elif len(bcp) == 0:
            solver_p_params = self.params.solver_p_neumann
        else:
            solver_p_params = self.params.solver_p_dirichlet

        for bc in bcp:
            bc.apply(A_p_corr)

        solver_p_corr = LinearSolver(*solver_p_params)
        solver_p_corr.set_operator(A_p_corr)
        if 'preconditioner' in solver_p_corr.parameters:
                solver_p_corr.parameters['preconditioner']['structure'] = 'same'
        solver_p_corr.parameters.update(self.params.p_corr_solver_parameters)

        timer.completed("create pressure correction solver")

        # Velocity correction solver
        if self.params.solver_u_corr not in ["WeightedGradient"]:
            # Velocity correction. Like for the tentative velocity, a single LHS is used.
            
            M = assemble(inner(u,v)*dx())
            A_u_corr = M
            if not self.params.low_memory_version and self.params.store_rhs_matrix_u_corr:
                rhs_u_corr = [None]*len(dims)
                Kp = [None]*len(dims)
                for d in dims:
                    Kp[d] = assemble(-k*inner(v, grad(p)[d])*dx())
                    rhs_u_corr[d] = RhsGenerator(U)
                    rhs_u_corr[d] += M, u1[d]
                    rhs_u_corr[d] += Kp[d], p1
                    rhs_u_corr[d] -= Kp[d], p0
            else:
                rhs_u_corr = [lambda: A_u_corr*u1[d].vector()+assemble(-k*inner(v, grad(p1-p0)[d])*dx()) for d in dims]

            # Apply BCs to LHS
            for bc in bcu:
                bc[0].apply(A_u_corr)

            solver_u_corr = LinearSolver(*self.params.solver_u_corr)
            solver_u_corr.set_operator(A_u_corr)
            if 'preconditioner' in solver_u_corr.parameters:
                solver_u_corr.parameters['preconditioner']['structure'] = 'same'
            solver_u_corr.parameters.update(self.params.u_corr_solver_parameters)

        elif self.params.solver_u_corr == "WeightedGradient":
            assert self.params.p_degree == 1
            #M1, _ = delta_memory()
            dPdX = _get_weighted_gradient(mesh, dims, v,p)
            #M2, _ = delta_memory()

            #import gc; gc.collect()
        """
        size = 0
        for m in dPdX:
            #print m.size(0), m.size(1)
            #for r in range(m.size(0)):
            for r in range(m.local_range(0)[0], m.local_range(0)[1]):
                row = m.getrow(r)
                size += row[0].nbytes+row[1].nbytes
        print "dPdX size: ", MPI.sum(mpi_comm_world(), size)/(1024*1024), "MB"
        print "M2-M1: ", M2-M1
        """

        timer.completed("create velocity correction solver")

        # Yield initial conditions
        yield ParamDict(spaces=spaces, observations=observations, controls=controls,
                        t=float(t), timestep=start_timestep, u=u1, p=p1)
        timer.completed("initial postprocessor update")

        # Time loop
        for timestep in xrange(start_timestep+1,len(timesteps)):
            t.assign(timesteps[timestep])

            # Update various functions
            problem.update(spaces, u1, p1, t, timestep, bcs, observations, controls)
            timer.completed("problem update")

            # Scale to solver pressure
            p0.vector()[:] *= 1.0/rho
            p1.vector()[:] *= 1.0/rho

            # Assemble convection
            if not self.params.low_memory_version and self.params.store_stiffness_matrix:
                # Assemble only convection matrix
                K_conv.zero()
                assemble(a_conv, tensor=K_conv)
                A_u_tent.axpy(-(1.0-theta), K_conv, True)
            else:
                # Assemble convection and diffusion matrix in one
                assemble(a2+a_conv, tensor=A)
            
            timer.completed("assemble convection matrix")
            
            # Build rhs for tentative velocity
            A_u_tent.axpy(-(1.0-theta), A, True)
            timer.completed("built A_u_tent for rhs")

            # Use A_u_tent in current form to create rhs
            # Note: No need to apply bcs to A_u_tent (this is set directly on b)
            b = [None]*len(dims)
            for d in dims:
                b[d] = rhs_u_tent[d]()
                for bc in bcu:
                    bc[d].apply(b[d])
            timer.completed("built tentative velocity rhs")

            # Construct lhs for tentative velocity
            A_u_tent.axpy(1.0, A, True)
            if not self.params.low_memory_version and self.params.store_stiffness_matrix:
                A_u_tent.axpy(1.0, K_conv, True)
            for bc in bcu:
                bc[0].apply(A_u_tent)

            timer.completed("u_tent construct lhs")

            # Check if preconditioner is to be rebuilt
            if timestep % self.params.rebuild_prec_frequency == 0 and 'preconditioner' in solver_u_tent.parameters:
                solver_u_tent.parameters['preconditioner']['structure'] = self.params.u_tent_prec_structure

            # Compute tentative velocity step
            for d in dims:
                iter = solver_u_tent.solve(A_u_tent, u1[d].vector(), b[d])

                # Preconditioner is the same for all three components, so don't rebuild several times
                if 'preconditioner' in solver_u_tent.parameters:
                    solver_u_tent.parameters['preconditioner']['structure'] = "same"

                timer.completed("u_tent solve (%s, %d dofs)"%(', '.join(self.params.solver_u_tent), b[d].size()), {"iter": iter})

            # Reset A_u_tent to mass matrix           
            A_u_tent.axpy(-theta, A, True)
            if not self.params.low_memory_version and self.params.store_stiffness_matrix:
                A_u_tent.axpy(-theta, K_conv, True)

            # Pressure correction
            b = rhs_p_corr()
            
            # Normalize around zero if no bcs are set on pressure
            # FIXME: Should really set nullspace
            if len(bcp) == 0:
                normalize(b)

            for bc in bcp:
                # Restore physical pressure and apply bcs
                b *= rho
                bc.apply(b)

                # Rescale to solver pressure
                b *= 1.0/rho

            timer.completed("p_corr construct rhs")

            # Solve p_corr
            iter = solver_p_corr.solve(p1.vector(), b)

            # Normalize around zero if no bcs are set on pressure
            # FIXME: Should really set nullspace
            if len(bcp) == 0:
                normalize(p1.vector())

            timer.completed("p_corr solve (%s, %d dofs)"%(', '.join(solver_p_params), b.size()), {"iter": iter})
            
            # Velocity correction
            if self.params.solver_u_corr not in ["WeightedGradient"]:
                for d in dims:
                    b = rhs_u_corr[d]()
                    for bc in bcu: bc[d].apply(b)
                    timer.completed("u_corr construct rhs")

                    iter = solver_u_corr.solve(u1[d].vector(), b)
                    timer.completed("u_corr solve (%s, %d dofs)"%(', '.join(self.params.solver_u_corr), b.size()),{"iter": iter})
            elif self.params.solver_u_corr == "WeightedGradient":
                for d in dims:
                    u1[d].vector().axpy(-dt, dPdX[d]*(p1.vector()-p0.vector()))
                    for bc in bcu:
                        bc[d].apply(u1[d].vector())
                    timer.completed("u_corr solve (weighted_gradient, %d dofs)" % u1[d].vector().size())

            # Update Adams-Bashford term for next timestep
            for d in dims:
                u_ab[d].vector().zero()
                u_ab[d].vector().axpy(1.5, u1[d].vector())
                u_ab[d].vector().axpy(-0.5, u0[d].vector())

             # Rotate functions for next timestep
            for d in dims:
                u0[d].assign(u1[d])
            p0.assign(p1)

            p0.vector()[:] *= rho
            p1.vector()[:] *= rho

            # Yield data for postprocessing
            yield ParamDict(spaces=spaces, observations=observations, controls=controls,
                            t=float(t), timestep=timestep, u=u1, p=p1, state=(u1,p1))
            timer.completed("updated postprocessing (completed timestep)")
