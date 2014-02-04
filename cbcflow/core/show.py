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

import dolfin
from .spaces import NSSpacePoolSplit
from cbcflow.dol import *

def animate_expression(f, name, V, t, timesteps):
    f0 = Function(V)
    for tv in timesteps:
        t.assign(tv)
        f0.interpolate(f)
        plot(f0, name="%s at t=%g" % (name, tv))

def animate_functions(functions, name, V):
    z = Function(V)
    for tk,zk in observations:
        z.assign(zk)
        plot(z, title="%s at t=%g" % (name, tk))

# TODO: Make version of this to test problem non-interactively, e.g. checking types and maybe some norms
def show_problem(problem, interactive=True, bc_snapshots=4):
    """Display properties of the problem.

    Intended for inspecting and debugging the problem setup.
    This functions runs through most of the interface
    """

    # Print params
    print("Problem parameters for problem of class %s:" % problem.__class__.__name__)
    print(str(problem.params))

    # Show the mesh
    plot(problem.mesh, title="Mesh")

    # Show eventual boundary markers from mesh
    if problem.facet_domains is not None:
        # TODO: Change max value to something just above the other values
        plot(problem.facet_domains, title="Facet domains")

    # Make time constant and animation timesteps
    t = Constant(problem.params.T0)
    n = bc_snapshots
    timesteps = [problem.params.T0 + (problem.params.T-problem.params.T0)*i/(n-1)
                 for i in range(n)]

    # Make linear function spaces suitable for plotting
    spaces = NSSpacePoolSplit(problem.mesh, u_degree=1, p_degree=1)
    V = spaces.V
    Q = spaces.Q

    # Plot observations if any
    observations = problem.observations(spaces, t)
    if isinstance(observations, list):
        if observations:
            animate_functions(observations, "Observation", observations[0][1].function_space(), t, timesteps)
    else:
        animate_expression(observations, "Observation", V, t, timesteps)

    # Plot controls if any
    controls = problem.controls(spaces)
    for i, c in enumerate(controls):
        plot(c, title="Control #%d" % (i,), mesh=problem.mesh)

    # Plot body force
    f = problem.body_force(spaces, t)
    f = as_vector(f)
    # TODO: Plot if not just zero
    #animate_expression(f, "Body force", V, t, timesteps)

    # Plot initial conditions
    u0, p0 = problem.initial_conditions(spaces, controls)
    u0 = as_vector(u0)
    plot(u0, title="Initial velocity", mesh=problem.mesh)
    plot(p0, title="Initial pressure", mesh=problem.mesh)

    # Get raw boundary conditions
    bcu, bcp = problem.boundary_conditions(spaces, u0, p0, t, controls)

    # Plot velocity BCs
    for i, tv in enumerate(timesteps):
        #problem.update(FIXME)
        for j, bc in enumerate(bcu):
            u, D = bc
            D = ("boundary domain %d " % D) if isinstance(D, int) else ""
            title = "Velocity BC #%d on %sat t=%g" % (j, D, t)
            u = as_vector(u)
            plot(u, title=title, mesh=problem.mesh)

    # Plot pressure BCs
    for i, tv in enumerate(timesteps):
        #problem.update(FIXME)
        for j, bc in enumerate(bcp):
            p, D = bc
            D = ("boundary domain %d " % D) if isinstance(D, int) else ""
            title = "Pressure BC #%d on %sat t=%g" % (j, D, t)
            plot(p, title=title, mesh=problem.mesh)

    if interactive:
        dolfin.interactive()
