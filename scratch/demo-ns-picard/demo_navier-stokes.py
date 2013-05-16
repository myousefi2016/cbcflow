"""This demo program solves the incompressible Navier-Stokes equations
on an L-shaped domain using Chorin's splitting method."""

# Copyright (C) 2010-2011 Anders Logg
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Mikael Mortensen 2011
#
# First added:  2010-08-30
# Last changed: 2011-06-30

# Begin demo

from dolfin import *
#from dolfin_adjoint import *


# Disable plotting:
if 0:
    def interactive(): pass
    def plot(*args,**kwargs): pass


# Print log messages only from the root process in parallel
parameters["std_out_all_processes"] = False;

# Load mesh from file
mesh = Mesh("lshape.xml.gz")
d = 2

class Inlet(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] > 1.0 - DOLFIN_EPS
class Outlet(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] > 1.0 - DOLFIN_EPS

boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(3)
DomainBoundary().mark(boundaries, 0)
Inlet().mark(boundaries, 1)
Outlet().mark(boundaries, 2)
ds = ds[boundaries]

n = FacetNormal(mesh)

# Define function spaces (P2-P1)
V1 = VectorFunctionSpace(mesh, "Lagrange", 1)
V2 = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)
M = V2*Q

# Define trial and test functions
u, p = TrialFunctions(M)
v, q = TestFunctions(M)

# Set parameter values
dt = 0.01
T = 3
nu = 0.01

# Define no-slip boundary conditions
bcs = [DirichletBC(M.sub(0), (0,)*d, boundaries, 0)]

# Define time-dependent pressure boundary condition
p_in = Expression("sin(3.0*t)", t=0.0)
t = p_in*n

# Create functions
up0 = Function(M) # Previous timestep
up1 = Function(M) # Last iterate of fixed-point in current timestep
u0, p0 = split(up0)
u1, p1 = split(up1)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

# Traction
#t = p*n - dot(n, nu*(grad(u) + grad(u).T)) # From equations
#t = Function(V1) # Make it a control living on boundary?
t_in = p_in*n # Make p_in a time-dependent control?
t_out = 0*n

# Picard linearization of Navier-Stokes, F = a*u - L = 0
a = k*(
    dot((1.0/k)*u + (grad(u)*u1), v)*dx

    + nu*inner(grad(u), grad(v))*dx
    #+ nu*dot(dot(n, grad(u).T), v)*ds()

    - p*div(v)*dx
    - q*div(u)*dx
    )
L = k*(
    dot((1.0/k)*u0 + f, v)*dx
    - dot(t_in, v)*ds(1)
    - dot(t_out, v)*ds(2)
    )
F = action(a, up1) - L

# Set initial condition
up0.interpolate(Expression(("0.0",)*(d+1)))
# Use initial condition as initial guess
up1.assign(up0)

# Use amg preconditioner if available
prec = "amg"

# Create solver
problem = NonlinearVariationalProblem(F, up1, bcs, J=a)
solver = NonlinearVariationalSolver(problem)

# Create files for storing solution
ufile = File("results/velocity.pvd")
pfile = File("results/pressure.pvd")

# Save initial condition to file
us, ps = up1.split()
ufile << us
pfile << ps

# Time-stepping
t = dt
while t < T + DOLFIN_EPS:

    # Update pressure boundary condition
    p_in.t = t

    # Solve for up1
    solver.solve()

    # Plot solution
    plot(us, title="Velocity", rescale=True)
    plot(ps, title="Pressure", rescale=True)

    # Save to file
    ufile << us
    pfile << ps

    # Move to next time step
    up0.assign(up1)
    t += dt
    print "t =", t

# Hold plot
interactive()

#assert replay_dolfin()
