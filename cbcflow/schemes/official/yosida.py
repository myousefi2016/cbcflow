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
from __future__ import division

__author__ = "Oyvind Evju <oyvinev@simula.no>"
__date__ = "2013-04-30"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from cbcflow.core.nsscheme import *
from cbcflow.core.rhsgenerator import *
from cbcflow.core.utils import Timer, epsilon, sigma, is_periodic
from cbcflow.core.timesteps import compute_regular_timesteps
from cbcflow.core.schemeutils import assign_ics_split, make_velocity_bcs, make_pressure_bcs, make_rhs_pressure_bcs
from cbcflow.core.spaces import NSSpacePoolSplit

from time import time

class Yosida(NSScheme):
    "Yosida scheme with lumped mass matrix in the Schur complement"

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.replace(
            # Default to P2-P1 (Taylor-Hood)
            u_degree = 2,
            p_degree = 1,
            )
        params.update(
            theta=0.5, # 0.5: Crank-Nicholson, 1.0: Backward Euler, 0.0: Forward Euler

            schur_tolerance=1e-8,
            conv_diff_tolerance=1e-8,
            supg = 0.0,
            )
        return params

    def solve(self, problem, update):
        from block import block_assemble
        from block.iterative import LGMRES
        from block.algebraic.trilinos.IFPACK import DD_ILU, DD_Jacobi, DD_Amesos
        from block.algebraic.trilinos.Epetra import LumpedInvDiag
        from block.algebraic.trilinos.MLPrec import ML

        # Get problem parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n  = FacetNormal(mesh)
        h = CellSize(mesh)

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Time(t0=timesteps[0])

        # Define function spaces
        spaces = NSSpacePoolSplit(mesh, self.params.u_degree, self.params.p_degree)
        V = spaces.V
        Q = spaces.Q

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # Functions
        u0 = Function(V, name="u0")
        u1 = Function(V, name="u1")
        u_corr = Function(V, name="u_corr")
        p0 = Function(Q, name="p0")
        p_corr = Function(Q, name="p_corr")
        
        # Get functions for data assimilation
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)

        # Get initial conditions
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_split(u0, p0, spaces, ics)
        u1.assign(u0)
        # Make scheme-specific representation of bcs
        bcs = problem.boundary_conditions(spaces, u0, p0, t, controls)
        bcu = make_velocity_bcs(problem, spaces, bcs)
        bcp = make_pressure_bcs(problem, spaces, bcs)
        Lbc = make_rhs_pressure_bcs(problem, spaces, bcs, v)

        # Problem coefficients
        nu = Constant(problem.params.mu/problem.params.rho)
        k  = Constant(dt)
        f  = as_vector(problem.body_force(spaces, t))
        theta = self.params.theta

        # Forms
        mu = inner(u,v)*dx()       
        au = inner(grad(u), grad(v))*dx()               
        u_conv = 0.5*(3*u1-u0)
        #u_conv = u1
        cu = inner(grad(u)*u_conv, v)*dx()
                            
        b1 = -inner(p, div(v))*dx()
        b2 = -inner(q, div(u))*dx()
        
        fu = 1/k*mu+theta*nu*au+theta*cu
        
        delta = Constant(self.params.supg)
        gamma = delta*dt*h
        supg = gamma*inner(grad(v)*u_conv, grad(u)*u_conv)*dx()
        
        Cu = assemble(cu, bcs=bcu)
        
        if float(delta) > 0:
            SUPG = assemble(supg, bcs=bcu)

        AA = block_assemble([[fu + supg, b1],[b2, 0]], bcs=[bcu, []])

        A = AA[0,0]
        D = AA[1,0]
        G = AA[0,1]
        
        L0 = ( (1/k)*dot(u1,v)*dx()
              - (1-theta)*dot(grad(u1)*u_conv, v)*dx()
              - (1-theta)*nu*inner(grad(u1), grad(v))*dx()
              #- (1-theta)*2*nu*inner(epsilon(u1), epsilon(v))*dx()
              #+ (1-theta)*nu*inner(grad(u1).T*n, v)*ds()
              #+ (1-theta)*p0*div(v)*dx()
              + Lbc
              + dot(f,v)*dx() )
        
        L1 = Constant(0)*q*div(u1)*dx()
        
        bb0 = assemble(L0, bcs=bcu)
        bb1 = assemble(L1)
        
        # Schur complement
        Mu = assemble(1/dt*mu, bcs=bcu)
        #Mu_inv = LumpedInvDiag(Mu)
        #Mu_inv = DD_ILU(Mu)
        #H = Mu_inv # First order (Chorin)
        #H = dt*Mu_inv - dt**2*Mu_inv*Mu_inv*Ku - dt**2*Mu_inv*Mu_inv*Cu # 2nd order, not always definite positive
        #H = dt*Mu_inv - dt**2*Mu_inv*(Ku+Cu) + dt**3*Mu_inv*(Mu_inv*(Ku+Cu))# 3rd order, always definite positive

        A_inv = DD_ILU(A)
        #A_inv = InvDiag(A)
        #A_inv = H
        #A_inv = DD_Amesos(A)
        
        S = D*A_inv*G # Schur complement
        
        # Preconditioner   
        Ap = assemble(dot(grad(p), grad(q))*dx())      
        
        for bc in bcp:
            bc.apply(Ap)

        Ap_inv = ML(Ap)
        X_inv = Ap_inv

        # Call update() with initial conditions
        update(u1, p0, float(t), start_timestep, spaces)

        guess = None
        # Loop over fixed timesteps
        for timestep in xrange(start_timestep+1,len(timesteps)):
            assign_time(t, timesteps[timestep])
        
            t0 = time()
            # Update various functions
            problem.update(spaces, u1, p0, t, timestep, bcs, observations, controls)

            t1 = time()
            print "Time spent updating: ", t1-t0
                
            # Update terms dependent on previous solution (convection terms)
            A.axpy(-1.0*theta, Cu, True)
            assemble(cu, bcs=bcu, tensor=Cu)
            A.axpy(1.0*theta, Cu, True)
            
            if float(delta) > 0:
                A.axpy(-1.0, SUPG, True)
                assemble(supg, bcs=bcu, tensor=SUPG)
                A.axpy(1.0, SUPG, True)
            
            assemble(L0, bcs=bcu, tensor=bb0)
            assemble(L1, tensor=bb1)
            
            t2 = time()
            print "Time spent re-assembling: ", t2-t1
            
            Aprec = DD_ILU(A)
            A_inv1 = LGMRES(A, precond=Aprec, initial_guess=u1.vector(), maxiter=250, nonconvergence_is_fatal=True, tolerance=self.params.conv_diff_tolerance)
            u1.vector()[:] = A_inv1*(bb0 - G*p0.vector())
            
            #S = D*Aprec*G # Schur complement (approximated)
            S_inv = LGMRES(S, precond=X_inv, tolerance=self.params.schur_tolerance, maxiter=100, show=1, nonconvergence_is_fatal=True)
            p_corr.vector()[:] = -S_inv * (bb1 - D * u1.vector())

            A_inv2 = LGMRES(A, precond=Aprec, initial_guess=u_corr.vector(), maxiter=250, nonconvergence_is_fatal=True, tolerance=self.params.conv_diff_tolerance)
            u_corr.vector()[:] = -A_inv2*G*p_corr.vector()           
            
            u1.vector()[:] += u_corr.vector()
            p0.vector()[:] += p_corr.vector()


            # Update last timestep
            u0.assign(u1)
            
            update(u1, p0, float(t), timestep, spaces)

        # Make sure annotation gets that the timeloop is over
        finalize_time(t)

        # Return some quantities from the local namespace
        states = (u1, p0)
        namespace = {
            "spaces": spaces,
            "observations": observations,
            "controls": controls,
            "states": states,
            "t": t,
            "timesteps": timesteps,
            }
        return namespace
