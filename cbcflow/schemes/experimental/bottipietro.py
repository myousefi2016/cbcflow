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

__author__ = "Martin Sandve Alnaes <kvs@simula.no>"
__date__ = "2013-08-13"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from ...core.nsscheme import *
from ...core.utils import Timer, is_periodic
from ...core.timesteps import compute_regular_timesteps
from ...core.schemeutils import (assign_ics_split,
                                make_velocity_bcs,
                                make_pressure_bcs,
                                make_penalty_pressure_bcs)
from ...core.spaces import NSSpacePoolSplit


class BottiPietro(NSScheme):
    "A pressure-correction scheme with discontinuous velocity and continuous pressure."

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.update(
            # Default to P2-P1
            u_degree = 2,
            p_degree = 1,
            )
        return params

    def solve(self, problem, update):
        # Notes about the original code:
        # - used CG3 for f
        # - had udeg=1 by default
        # - had a hack for facetarea

        # Get problem parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n  = FacetNormal(mesh)
        cell = mesh.ufl_cell()
        x = cell.x

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Time(t0=timesteps[start_timestep])

        # Define function spaces
        spaces = NSSpacePoolSplit(mesh, self.params.u_degree, self.params.p_degree, u_family="DG")
        V = spaces.V
        Q = spaces.Q
        udeg = self.params.u_degree

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
        ubcs, _pbcs = bcs
        ubcs = [(as_vector(function), region) for function, region in ubcs]
        #bcu = make_velocity_bcs(problem, spaces, bcs)
        bcp = make_pressure_bcs(problem, spaces, bcs)

        # Remove boundary stress term if problem is periodic
        beta = 0 if is_periodic(bcp) else 1

        # Problem coefficients
        nu = Constant(problem.params.mu/problem.params.rho, name="nu")
        k  = Constant(dt, name="dt")
        kinv = 1.0 / k
        f  = as_vector(problem.body_force(spaces, t))



        # Define penalty terms
        eta = Constant(5.0, cell=cell)
        np = n('+')
        h_T = CellSize(mesh)
        h_dT = FacetArea(mesh)

        h_F_ext = h_T / h_dT
        if 1:
            h_T_min = Min(h_T('+'), h_T('-'))
            h_F_int = h_T_min / h_dT('+') # TODO: Don't require restricting this in UFL/FFC, not necessary.
        else:
            h_F_int = avg(h_T)

        # Penalty term
        pen_ext = eta * udeg**2 / h_F_ext
        pen_int = eta('+') * udeg**2 / h_F_int

        # Free indices for use with implicit summation
        i, j = indices(2)

        # Left hand side for u equations
        a_u  = inner(grad(u), grad(v))*dx()
        a_u += ( pen_int*dot(jump(u),jump(v))
                - dot( avg(grad(u))*np, jump(v) )
                - dot( avg(grad(v))*np, jump(u) ) ) * dS()
        a_u += ( pen_ext*dot(u,v)
                - dot(grad(u)*n, v)
                - dot(grad(v)*n, u) ) * ds()
        a_u += kinv*dot(u, v)*dx() # Time derivative term

        # Right hand side for u equations
        b_u = -dot(f,v)*dx()                       # Forcing term
        b_u += kinv*dot(u0,v)*dx()                 # Time derivative term
        b_u -= dot(v, 2*grad(p1) - grad(p0))*dx() # Pressure coupling term

        for ubc, region in ubcs:
            if 0: # FIXME: What is the right formulation for this term?
                b_u += pen_ext*dot(ubc,v)*ds(region)
            else:
                b_u += pen_ext*dot(ubc,n)*dot(v,n)*ds(region)
            b_u -= dot(grad(v)*n, ubc)*ds(region)

        # Pressure Poisson equation
        a_p = dot(grad(p),grad(q))*dx()
        b_p  = dot(grad(p0), grad(q))*dx()
        b_p -= kinv*div(u1)*q*dx()
        b_p -= kinv('+')*dot(np, jump(u1))*avg(q)*dS()
        for ubc, region in ubcs:
            b_p += kinv*dot(n, u1-ubc)*q*ds(region)

        # Assemble time independent matrices
        A_u = assemble(a_u)
        A_p = assemble(a_p)

        # Apply BCs to matrices
        for bc in bcp:
            bc.apply(A_p)

        if 0:
            solver_name = "lu"
            solver_params = None
        else:
            solver_name = "gmres"
            solver_params = {
                'relative_tolerance': 1e-15,
                'monitor_convergence': True,
                'gmres': { 'restart': 100 },
                }

        solver_u = LinearSolver("bicgstab", "hypre_euclid")
        solver_u.set_operator(A_u)
        #solver_u.parameters['preconditioner']['reuse'] = True

        solver_p = LinearSolver("gmres", "ml_amg")
        solver_p.set_operator(A_p)
        #solver_p.parameters['preconditioner']['reuse'] = True

        # Call update() with initial conditions
        update(u0, p0, float(t), 0, spaces)

        timer = Timer(self.params.enable_timer)

        # Loop over fixed timesteps
        for timestep in xrange(start_timestep+1,len(timesteps)):
            assign_time(t, timesteps[timestep])

            # Update various functions
            problem.update(spaces, u0, p0, t, timestep, bcs, observations, controls)
            timer.completed("problem update")

            # Solve advection-diffusion equations
            bh_u = assemble(b_u)
            solver_u.solve(u1.vector(), bh_u)

            # Solve pressure equation
            p0.assign(p1)
            bh_p = assemble(b_p)
            #for bc in bcp: bc.apply(A_p, bh_p)
            for bc in bcp:
                bc.apply(bh_p)
            solver_p.solve(p1.vector(), bh_p)

            # Compute change in u from last timestep
            #uchange = assemble((u1-u0)**2*dx())
            #print "uchange:", uchange

            # Update postprocessing
            update(u1, p1, float(t), timestep, spaces)

            # Rotate functions for next timestep
            u0.assign(u1)
            #p0.assign(p1) # Must wait until the next advection-diffusion step is done


        # Make sure annotation gets that the timeloop is over
        finalize_time(t)

        # Return some quantities from the local namespace
        states = (u1, p1)
        namespace = {
            "spaces": spaces,
            "observations": observations,
            "controls": controls,
            "states": states,
            "t": t,
            "timesteps": timesteps,
            }
        return namespace
