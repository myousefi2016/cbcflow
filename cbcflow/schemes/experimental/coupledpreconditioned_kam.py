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
from cbcflow.utils.common import Timer, epsilon, sigma, is_periodic
from cbcflow.utils.schemes import (RhsGenerator,
                                   compute_regular_timesteps,
                                   assign_ics_split,
                                   make_velocity_bcs,
                                   make_pressure_bcs,
                                   make_rhs_pressure_bcs)
from cbcflow.utils.core import NSSpacePoolSplit

from time import time


class CoupledPreconditonedKAM(NSScheme):
    "Coupled scheme with block preconditioning using cbc.block"

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
            theta=1.0,#0.5, # 0.5: Crank-Nicholson, 1.0: Backward Euler, 0.0: Forward Euler

            fixed_point_tolerance=1e-6,
            max_fixed_point_iterations=500,
            )
        return params

    def solve(self, problem, update):
        from block import block_assemble, block_vec, block_mat
        from block.iterative import LGMRES
        from block.algebraic.trilinos.IFPACK import DD_ILU, DD_Jacobi
        from block.algebraic.trilinos.Epetra import LumpedInvDiag
        from block.algebraic.trilinos.MLPrec import ML

        # Get problem parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n  = FacetNormal(mesh)

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Time(t0=timesteps[start_timestep])

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
        p0 = Function(Q, name="p0")
        p1 = Function(Q, name="p1")

        # Get functions for data assimilation
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)

        # Get initial conditions
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_split(u0, p0, spaces, ics)
        u1.assign(u0)
        p1.assign(p0)

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
        mu0 = inner(u0,v)*dx()

        mp = p*q*dx()
        mp0 = p0*q*dx()
        mpmuinv = (1.0/nu)*p*q*dx()

        #au = 2*inner(epsilon(u), epsilon(v))*dx()
        #au0 = 2*inner(epsilon(u0), epsilon(v))*dx()

        au = inner(grad(u), grad(v))*dx()
        au0 = inner(grad(u0), grad(v))*dx()


        ap = inner(grad(p), grad(q))*dx()
        ap0 = inner(grad(p0), grad(q))*dx()

        cu = inner(grad(u)*u0, v)*dx()
        cu0 = inner(grad(u0)*u0, v)*dx()

        cp = inner(grad(p),u0)*q*dx()
        cp0 = inner(grad(p0), u0)*dx()


        b1 = -inner(p, div(v))*dx()
        b10 = -inner(p0, div(v))*dx()

        b2 = -inner(q, div(u))*dx()
        b20 = -inner(q, div(u0))*dx()

        fu = 1/k*mu+theta*nu*au+theta*cu
        fp = 1/k*mp+theta*nu*ap+theta*cp


        Cu = assemble(cu, bcs=bcu)

        # Preconditioner
        Lp = Constant(0)*q*dx()
        Ap, _ = assemble_system(k*ap, Lp, bcp)
        Fp = assemble(fp)
        Mp = assemble(mp)
        Mpmuinv = assemble(mpmuinv)
        Fu = assemble(fu)
        B1 = assemble(b1)
        B2 = assemble(b2)
        Cp = assemble(k*cp)
        #Cu = assemble(cu)


        for bc in bcu:
            bc.apply(Fu)
            #bc.apply(B1)

#        for bc in bcp:
#            bc.apply(Ap)
#            bc.apply(Fp)
#            bc.apply(Mp)


        AA = block_assemble([[fu, theta*b1],[theta*b2, 0]], bcs=[bcu, []])


        L0 = ( (1/k)*dot(u0,v)*dx()
              - (1-theta)*dot(grad(u0)*u0, v)*dx()
              - (1-theta)*nu*inner(grad(u0), grad(v))*dx()
              + (1-theta)*p0*div(v)*dx()
              + (1-theta)*q*div(u0)*dx()
              + Lbc
              + dot(f,v)*dx() )
        L1 = Constant((1-theta))*q*div(u0)*dx()


        bb0 = assemble(L0, bcs=bcu)
        bb1 = assemble(L1)

        bb = block_vec([bb0,bb1])


        Fp_inv = DD_ILU(Fp)
#        Mp_inv = DD_Jacobi(Mp)
        Mp_inv = LumpedInvDiag(Mp)
        Mpmuinv_inv = LumpedInvDiag(Mpmuinv)
        Ap_inv = ML(Ap)

        #Fu_inv = DD_Amesos(Fu)
        #Fp_inv = DD_Amesos(Fp)
        #Mp_inv = DD_Amesos(Mp)
        #Ap_inv = DD_Amesos(Ap)


        X = Ap*Fp_inv*Mp
        X_inv = Ap_inv
        X_inv2 = Mpmuinv_inv

#        P_inv = block_mat([[Fu_inv, 0], [0, 1]])
#        P_inv *= block_mat([[1, -B1], [0, 1]])
#        P_inv *= block_mat([[1, 0],[0, X_inv]])
#        P_inv = block_mat([[Fu_inv, 0], [0, X_inv]])


#        P_inv = P_inv.block_collapse()
#        P_inv.scheme('sgs')

        '''

        P_inv = block_mat([[Fu_inv, 0],[0, 1]])#.scheme('sgs')
        P_inv *= block_mat([[1, -B1], [0, 1]])#.scheme('sgs')
        P_inv *= block_mat([[1, 0], [0, -Mp_inv]])#.scheme('sgs')
        P_inv *= block_mat([[1, 0], [0, Fp]])#.scheme('sgs')
        P_inv *= block_mat([[1, 0], [0, Ap_inv]])#.scheme('sgs')
        '''
        '''
        P_inv = block_mat([[Fu_inv, 0],[0, 1]]).#scheme('sgs')
        P_inv *= block_mat([[1, -B1], [0, 1]])#.scheme('sgs')
        P_inv *= block_mat([[1, 0], [0, -Mp_inv]])#.scheme('sgs')
        P_inv *= block_mat([[1, 0], [0, Fp]])#.scheme('sgs')
        P_inv *= block_mat([[1, 0], [0, Ap_inv]])#.scheme('sgs')
        '''

        guess = block_vec([u0.vector().copy(), p0.vector().copy()])

        # Call update() with initial conditions
        update(u0, p0, float(t), start_timestep, spaces)

        # Loop over fixed timesteps
        for timestep in xrange(start_timestep+1,len(timesteps)):
            assign_time(t, timesteps[timestep])

            # Update various functions
            problem.update(spaces, u0, p0, t, timestep, bcs, observations, controls)
            t1 = time()

            print norm(u0)

            # Update terms dependent on previous solution (convection terms)
            AA[0,0].axpy(-1.0*theta, Cu, True)
            Fu.axpy(-1.0*theta, Cu, True)
            assemble(cu, bcs=bcu, tensor=Cu)
            AA[0,0].axpy(1.0*theta, Cu, True)
            Fu.axpy(1.0*theta, Cu, True)

            Fu_inv = DD_ILU(Fu)
            P_inv = block_mat([[Fu_inv, 0], [0, X_inv + X_inv2]])


            Fp.axpy(-1.0, Cp, True)
            assemble(cp, tensor=Cp)
            Fp.axpy(1.0, Cp, True)

            assemble(L0, bcs=bcu, tensor=bb0)
            assemble(L1, tensor=bb1)

            bb[0] = bb0
            bb[1] = bb1

            guess[0] = u0.vector().copy()
            guess[1] = p0.vector().copy()

            #guess = block_vec([u0.vector().copy(), p0.vector().copy()])

            t2 = time()

            print "Time spent assembling: ", t2-t1

            t1 = time()
            TOL = 1e-6
            AAinv = LGMRES(AA, precond=P_inv, tolerance=TOL, maxiter=500, show=1, initial_guess=guess)

            u_vector, p_vector = AAinv * bb
            t2 = time()
            print "Time spent solving: ", t2-t1
            u1.vector()[:] = u_vector
            p1.vector()[:] = p_vector

            # Update last timestep
            u0.assign(u1)
            p0.assign(p1)

            #print "Fixed point iteration converged in %d iterations. (err=%.4e)" %(count, err)
            update(u0, p0, float(t), timestep, spaces)

        # Make sure annotation gets that the timeloop is over
        finalize_time(t)

        # Return some quantities from the local namespace
        states = (u0, p0)
        namespace = {
            "spaces": spaces,
            "observations": observations,
            "controls": controls,
            "states": states,
            "t": t,
            "timesteps": timesteps,
            }
        return namespace
