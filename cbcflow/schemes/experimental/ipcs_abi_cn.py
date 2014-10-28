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


from cbcflow.core.nsscheme import *
from cbcflow.utils.common import is_periodic, epsilon
from cbcflow.utils.schemes import (RhsGenerator,
                                   compute_regular_timesteps,
                                   assign_ics_segregated,
                                   make_segregated_velocity_bcs,
                                   make_pressure_bcs,
                                   make_penalty_pressure_bcs)
from cbcflow.utils.core import NSSpacePoolSegregated

from dolfin import *
#from os import getpid, path, makedirs, getcwd
import os
from commands import getoutput
import time as tme
from numpy import ceil

#####################################################################
#from solverbase import *

#master=MPI.process_number()==0
#parameters["std_out_all_processes"] = False;

class IPCS_ABI_CN(NSScheme):
    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.update(
            # Default to P1-P1
            u_degree = 1,
            p_degree = 1,
            #theta = 0.5,
            #rebuild_prec_frequency = 1e16,
            #u_tent_prec_structure = "same_nonzero_pattern",
            #u_tent_solver_parameters = {},
            #p_corr_solver_parameters = {},
            #u_corr_solver_parameters = {},
            )
        return params

    def solve(self, problem, update, timer):
        # Get problem parameters
        mesh     = problem.mesh  #Mesh
        dt       = problem.params.dt    #Time step size
        t        = problem.params.T0     #Current time
        T        = problem.params.T     #End time
        tstep    = t/T           #Time steps
        max_iter = 1             #Pressure velocity
        max_error =.0001         # Percent of norm
        iter_1   = 1           #Pressure velocity on 1st TS

        #Solver tolerances
        at  = 1e-15   #Krylov solver absolute tolerance
        rt  = 1e-15   #Krylov solver relative tolerance
        pat = 1e-15   #Krylov solver abolute  tolerance for pressure
        prt = 1e-15   #Krylov solver relative tolerance for pressure

        #Some other stuff
        dim     = mesh.geometry().dim()               #Dimension of geometry
        f       = Constant((0,)*dim)                  #Body force
        nu      = Constant(problem.params.mu/problem.params.rho)                #kinametic Viscosity
        #element order for pressure
        #con_dom = self.options["constrained_domain"]  #constrained domain
        uOrder = self.params.u_degree
        pOrder = self.params.p_degree

        #Solver Parameters
        parameters["linear_algebra_backend"]        = "PETSc"
        parameters["form_compiler"]["optimize"]     = False
        parameters["form_compiler"]["cpp_optimize"] = True
        set_log_active(True)

        # Declare solution Functions and FunctionSpaces
        spaces = NSSpacePoolSegregated(mesh, self.params.u_degree, self.params.p_degree)

        V = spaces.U
        VV = spaces.V
        Q = spaces.Q
        #V = FunctionSpace      (mesh,'CG',uOrder)
        #VV= VectorFunctionSpace(mesh,'CG',uOrder)
        #Q = FunctionSpace      (mesh,'CG',pOrder)
        u = TrialFunction(V)
        v = TestFunction (V)
        p = TrialFunction(Q)
        q = TestFunction (Q)

        #Assign array for velocity and pressure
        if dim == 2: u_components = ['u0', 'u1']
        else:        u_components = ['u0', 'u1', 'u2']
        sys_comp  =  u_components + ['p']

        # Use dictionaries to hold all Functions and FunctionSpaces
        VV = dict((ui, V) for ui in u_components); VV['p'] = Q

        # Start from previous solution if restart_folder is given

        ics=problem.initial_conditions(spaces, None)

        #import ipdb; ipdb.set_trace()


        q_     = dict((ui, interpolate(ics[0][u_components.index(ui)],V)) for ui in u_components)
        q_1    = dict((ui, interpolate(ics[0][u_components.index(ui)],V)) for ui in u_components)
        q_2    = dict((ui, interpolate(ics[0][u_components.index(ui)],V)) for ui in u_components)

        q_['p']= interpolate(ics[-1],Q)


        u_  = as_vector([q_[ui]  for ui in u_components]) # Velocity vector at t
        u_1 = as_vector([q_1[ui] for ui in u_components]) # Velocity vector at t - dt
        u_2 = as_vector([q_2[ui] for ui in u_components]) # Velocity vector at t - 2*dt
        p_  = q_['p']                # pressure at t - dt/2
        dp_ = Function(Q)            # pressure correction

        #Degrees of Freedom
        num_dofs = sum(([q_[ui].vector().size()  for ui in u_components]))
        tot_num_dofs = sum(([q_[ui].vector().size()  for ui in sys_comp]))

        print '-'*60
        print 'Degrees of freedom for velocity are:          ',     num_dofs
        print 'Degrees of freedom for velocity and pressure: ', tot_num_dofs

#-----------------------------------------------------------------------------------------
#                             Krylov Solver Tolerances
#-----------------------------------------------------------------------------------------

        print '-'*60
        print 'Krylov solver relative tolerance is ', rt
        print 'Krylov solver absolute tolerance is ', at
        print 'Krylov solver relative tolerance for the pressure is ', prt
        print 'Krylov solver absolute tolerance for the pressure is ', pat
        print '-'*60

        u_sol = KrylovSolver('bicgstab', 'jacobi')
        u_sol.parameters['error_on_nonconvergence']  = False
        u_sol.parameters['nonzero_initial_guess']    = True
        u_sol.parameters['monitor_convergence']      = False
        u_sol.parameters['relative_tolerance']       = rt
        u_sol.parameters['absolute_tolerance']       = at
        reset_sparsity = True

        du_sol = KrylovSolver('bicgstab', 'jacobi')
        du_sol.parameters['error_on_nonconvergence'] = False
        du_sol.parameters['nonzero_initial_guess']   = True
        du_sol.parameters['preconditioner']['structure'] = "same"
        #du_sol.parameters['monitor_convergence']    = True
        du_sol.parameters['relative_tolerance']      = rt
        du_sol.parameters['absolute_tolerance']      = at

        p_sol = KrylovSolver('gmres', 'hypre_amg')
        p_sol.parameters['error_on_nonconvergence']  = False
        p_sol.parameters['nonzero_initial_guess']    = True
        p_sol.parameters['preconditioner']['structure']  = "same"
        #p_sol.parameters['monitor_convergence']     = True
        p_sol.parameters['relative_tolerance']       = prt
        p_sol.parameters['absolute_tolerance']       = pat


#-----------------------------------------------------------------------------------------
#                             Boundary Conditions
#-----------------------------------------------------------------------------------------

        # Get initial and boundary conditions
        #bcs = problem.boundary_conditions(V, Q, t)
        bcs = problem.boundary_conditions(spaces, u_, p_, t, None)

        #bcu, bcp = bcs[:-1], bcs[-1]
        #bcu, bcp = bcs
        bcu = make_segregated_velocity_bcs(problem, spaces, bcs)
        bcp = make_pressure_bcs(problem, spaces, bcs)
        bcs = dict((ui, []) for ui in sys_comp)
        #import ipdb; ipdb.set_trace()
        #Assign Boundary Conditions to bc array
        if dim == 2:
          #bcs['u0'] = bcu[0]
          bcs['u0'] = [_bc[0] for _bc in bcu]
          #bcs['u1'] = bcu[1]
          bcs['u1'] = [_bc[1] for _bc in bcu]
          bcs['p'] =  bcp[:]
        else:
          bcs['u0'] = bcu[0]
          bcs['u1'] = bcu[1]
          bcs['u2'] = bcu[2]
          bcs['p'] =  bcp[:]

#-----------------------------------------------------------------------------------------
#                             Defining Equations
#-----------------------------------------------------------------------------------------

        # Preassemble some constant in time matrices
        M = assemble(inner(u, v)*dx)                 # Mass matrix
        K = assemble(nu*inner(grad(u), grad(v))*dx)  # Diffusion matrix
        Ap = assemble(inner(grad(q), grad(p))*dx)    # Pressure Laplacian
        A = Matrix()                                 # Coefficient matrix (needs reassembling)

        # Apply boundary conditions on M and Ap that are used directly in solve
        [bc.apply(M)  for bc in bcs['u0']]
        [bc.apply(Ap) for bc in bcs['p']]

        #Compress Pressure Laplacian
        try:    Ap.compress()
        except: pass

        # Adams Bashforth projection of velocity at t - dt/2
        U_ = 1.5*u_1 - 0.5*u_2

        # Convection form
        a  = 0.5*inner(v, dot(U_, nabla_grad(u)))*dx

        # Preassemble constant pressure gradient matrix
        P = dict((ui, assemble(v*p.dx(i)*dx)) for i, ui in enumerate(u_components))

        # Preassemble velocity divergence matrix
        if uOrder==pOrder: R = P
        else: R = dict((ui, assemble(q*u.dx(i)*dx)) for i, ui in  enumerate(u_components))


        x_  = dict((ui, q_ [ui].vector()) for ui in sys_comp)     # Solution vectors t
        x_1 = dict((ui, q_1[ui].vector()) for ui in u_components) # Solution vectors t - dt
        x_2 = dict((ui, q_2[ui].vector()) for ui in u_components) # Solution vectors t - 2*dt
        b   = dict((ui, Vector(x_[ui])) for ui in sys_comp)       # rhs vectors
        bold= dict((ui, Vector(x_[ui])) for ui in sys_comp)       # rhs temp storage vectors
        work = Vector(x_['u0'])


#-----------------------------------------------------------------------------------------
#                             Time Loop
#-----------------------------------------------------------------------------------------

        #Time loop parameters
        dt_         = dt #Time step size
        total_iters = 0  #total iterations

        #Logs
        #self.start_timing()
        set_log_active(False)

        while t < T:
            t     += dt_      #update current time
            tstep += 1        #update time step
            j      = 0        #Iterations
            err = 1e8         #Initial error
            total_iters += 1  #update total_iterations

            problem.update(spaces, u_, p_, t, tstep, bcs, None, None)

            #Use iter_1 on 1st 10 time steps
            if tstep < 10:   num_iter = max(iter_1, max_iter)
            else:            num_iter = max_iter

            #Loop over inner iterations
            while err > max_error and j < num_iter:

                j += 1         #update number of iterations

                ### Start by solving for an intermediate velocity ###
                if j == 1:

                    # Only on the first iteration because nothing here is changing in time
                    # Set up coefficient matrix for computing the rhs:
                    # Warning! Must reset A = Matrix() for periodic boundary conditions(?)
                    A = assemble(a, tensor=A)
                    A._scale(-1.)            # Negative convection on the rhs
                    A.axpy(1./dt_, M, True)  # Add mass
                    A.axpy(-0.5, K, True)    # Add diffusion

                    # Compute rhs for all velocity components
                    for ui in u_components:
                        b[ui][:] = A*x_1[ui]

                    # Reset matrix for lhs
                    A._scale(-1.)
                    A.axpy(2./dt_, M, True)
                    [bc.apply(A) for bc in bcs['u0']]

                    print "Matrix (tentative) norm: ", A.norm('frobenius')
                    for i in range(2):
                        print (1.5*u_1[i].vector()[:]-0.5*u_2[i].vector()[:]).norm('l2')
                    #import ipdb; ipdb.set_trace()
                    #exit()

                    #Solving for tentative velocity

                    for ui in u_components:
                        bold[ui][:] = b[ui][:]
                        b[ui].axpy(-1., P[ui]*x_['p'])
                        [bc.apply(b[ui]) for bc in bcs[ui]]
                        #u_sol.solve(A, x_[ui], b[ui])
                        solve(A, x_[ui], b[ui])
                        print "Rhs (tentative) norm: ", b[ui].norm('l2')
                        print "ui norm: ", x_[ui].norm('l2')
                        b[ui][:] = bold[ui][:]


                ### Solve pressure ###
                dp_.vector()[:] = x_['p'][:]
                b['p'][:] = Ap*x_['p']
                for ui in u_components:
                  b['p'].axpy(-1./dt_, R[ui]*x_[ui]) # Divergence of u_
                [bc.apply(b['p']) for bc in bcs['p']]
                p_sol.solve(Ap, x_['p'], b['p'])
                if normalize: normalize(x_['p'])
                dp_.vector()[:] = x_['p'][:] - dp_.vector()[:]

                oldnorm = 0; newnorm = 0
                ### Update velocity ###
                for ui in u_components:
                    b[ui][:] = M*x_[ui][:]
                    b[ui].axpy(-dt_, P[ui]*dp_.vector())
                    [bc.apply(b[ui]) for bc in bcs[ui]]
                    oldnorm += norm(x_[ui][:])
                    du_sol.solve(M, x_[ui], b[ui])
                    newnorm += norm(x_[ui][:])
                try:
                    err = 100*abs(oldnorm - newnorm)/oldnorm
                except:
                        if master:
                                print

                #if master:
                #        print 'Percent change of norm for velocity update', err, 'Iter ',j
            # Update to a new timestep
            for ui in u_components:
                x_2[ui][:] = x_1[ui][:]
                x_1[ui][:] = x_ [ui][:]
            #if master:
            #    print
            ################ Hack!! Because PETSc bicgstab with jacobi errors on the first tstep and exits in parallel ##
            if tstep == 0:
                u_sol = KrylovSolver('bicgstab', 'jacobi')
                u_sol.parameters['error_on_nonconvergence'] = False
                u_sol.parameters['nonzero_initial_guess']   = True
                u_sol.parameters['monitor_convergence']     = True
    #################################################################################################


           # Update
            #plot(u_)
            #interactive()
            #import ipdb; ipdb.set_trace()
            for i, ui in enumerate(u_):
                print "u[%d] norm: " %i, ui.vector().norm('l2')
            print "p norm: ", p_.vector().norm('l2')

            update(u_, p_, float(t), tstep, spaces)
            #self.update_temp(problem, t, u_, q_['p'])


        # Return some quantities from the local namespace
        states = (u_, p_)
        namespace = {
            "spaces": spaces,
            "observations": None,
            "controls": None,
            "states": states,
            "t": t,
            "timesteps": tstep,
            }
        return namespace

    def __str__(self):
        name = "IPCS_opt_AB_CN"
        return name
