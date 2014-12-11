from dolfin import *

dt = 1e-1
nx = 100

mesh = UnitSquareMesh(nx, nx)
d = mesh.geometry().dim()
n = FacetNormal(mesh)

V = VectorFunctionSpace(mesh, "CG", 1)
U = FunctionSpace(mesh, "CG", 1)
Q = U #FunctionSpace(mesh, "CG", 1)

v = [TestFunction(U) for i in range(d)]
q = TestFunction(Q)
ut = [TrialFunction(U) for i in range(d)]
pt = TrialFunction(Q)

rho = Constant(1.06, name="rho")
mu = Constant(1.0/30.0, name="mu")
nu = Constant(float(mu)/float(rho), name="nu")

k = Constant(dt, name="k")
theta = Constant(0.5, name="theta")

# Setting alpha = 0 gives mass matrix in velocity correction (dropping stiffness term)
#alpha = 1
alpha = 0

def SegregatedFunction(U, dim, name):
    return as_vector([Function(U, name="{0}[comp={1}]".format(name, i)) for i in range(dim)])

def SegregatedTimeFunction(U, dim, name, timesteps):
    return { ts: SegregatedFunction(U, dim, "{0}[ts={1}]".format(name, ts)) for ts in timesteps }

def TimeFunction(Q, name, timesteps):
    return { ts: Function(Q, name="{0}[ts={1}]".format(name, ts)) for ts in timesteps }

# Input functions with history
f = SegregatedTimeFunction(U, d, "f", (0,+1))
g = SegregatedTimeFunction(U, d, "g", (0,+1))

# Solution functions with or without history
u = SegregatedTimeFunction(U, d, "u", (-1,0,+1))
utent = SegregatedFunction(U, d, "utent")
p = TimeFunction(Q, "p", (0, 1))

# Extrapolate from solutions at t[-1], t[0] to t=0.5*(t[0]+t[1])
u_ext = as_vector([1.5*u[0][i] + 0.5*u[-1][i] for i in range(d)])
#p_ext = 1.5*p[0] + 0.5*p[-1]
p_ext = p[0]

a_u_tent = [(
    (1.0/k) * inner(ut[i], v[i])*dx
    + theta * inner(dot(grad(ut[i]), u_ext), v[i])*dx
    + theta*nu * inner(grad(ut[i]), grad(v[i]))*dx # (iii)
    - theta*nu * inner(dot(grad(ut[i]), n), v[i])*ds # (iii)
    ) for i in range(d)]

b_u_tent = [(
    (1.0/k)*inner(u[0][i], v[i])*dx
    - (1-theta) * inner(dot(grad(u[0][i]), u_ext), v[i])*dx
    - (1-theta)*nu * inner(grad(u[0][i]), grad(v[i]))*dx # (ii)
    + (1.0/rho) * p_ext*v[i].dx(i)*dx # (i)
    + inner(f[0][i], v[i])*dx
    # Boundary terms from partial integration:
    - (1.0/rho) * inner(p_ext*n[i], v[i])*ds # (i)
    - (1-theta)*nu * inner(dot(grad(u[0][i]), n), v[i])*ds # (ii)
    ) for i in range(d)]

a_p = (
    (1.0/rho) * inner(grad(pt), grad(q))*dx # (i)
    # Boundary terms from partial integration:
    - (1.0/rho) * inner(dot(grad(pt), n), q)*ds # (i)
    )
b_p = (
    (1.0/rho) * inner(grad(p_ext), grad(q))*dx # (ii)
    #+ alpha*(theta*nu) * inner(grad(div(utent)), grad(q))*dx # (iii) for u degree 2, grad(div(utent)) is nonzero, can we drop it then?
    - (1.0/k)*div(utent)*q*dx
    #+ (1.0/k)*inner(utent, grad(q))*dx (iv)
    #- (1.0/k)*inner(dot(utent, n), q)*ds (iv)
    # Boundary terms from partial integration:
    - (1.0/rho) * inner(dot(grad(p_ext), n), q)*ds # (ii)
    #- alpha*(theta*nu) * inner(dot(grad(div(utent)), n), q)*ds # (iii) for u degree 2, grad(div(utent)) is nonzero, can we drop it then?
    )

a_u_corr = [(
    inner(ut[i], v[i])*dx
    + alpha*(k*theta*nu) * inner(grad(ut[i]), grad(v[i]))*dx # (i)
    # Boundary terms from partial integration:
    - alpha*(k*theta*nu) * inner(dot(grad(ut[i]), n), v[i])*ds # (i)
    ) for i in range(d)]
b_u_corr = [(
    inner(u[0][i], v[i])*dx
    + alpha*(k*theta*nu) * inner(grad(utent[i]), grad(v[i]))*dx # (ii)
    + k/rho * inner(p[+1].dx(i), v[i])*dx
    - k/rho * inner(p_ext.dx(i), v[i])*dx
    # Boundary terms from partial integration:
    + alpha*(k*theta*nu) * inner(dot(grad(utent[i]), n), v[i])*ds # (ii)
    ) for i in range(d)]

# Setup strong bcs
bcs_u_tent = []
bcs_p = []
bcs_u_corr = []

# Preassemble
lhs_p = assemble(a_p)
rhs_p = assemble(b_p)
for bc in bcs_p:
    bc.apply(lhs_p)

lhs_u_corr = [assemble(a_u_corr[i]) for i in range(d)]
rhs_u_corr = [assemble(b_u_corr[i]) for i in range(d)]
for i in range(d):
    for bc in bcs_u_corr:
        bc.apply(lhs_u_corr[i])

t = 0.0
T = dt
while t < T:
    # Set t0,t to the beginning and end of the
    # interval we're computing over now
    t0 = float(t)
    t += dt

    # Set u^-1 to be u at t0-dt
    for i in range(d):
        u[-1][i].assign(u[0][i])

    # Set u^0 to be u at t0
    for i in range(d):
        u[0][i].assign(u[1][i])

    # Find tentative u
    lhs_u_tent = [assemble(a_u_tent[i]) for i in range(d)]
    rhs_u_tent = [assemble(b_u_tent[i]) for i in range(d)]
    for i in range(d):
        for bc in bcs_u_tent:
            bc.apply(rhs_u_tent[i])
    for i in range(d):
        solve(lhs_u_tent[i], rhs_u_tent[i], utent[i].vector()) # TODO: Setup solvers outside of loop

    # Find p^1 (pressure at time t)
    for bc in bcs_p:
        bc.apply(rhs_p)
    solve(lhs_p, rhs_p, p[1].vector()) # TODO: Setup solvers outside of loop

    # Find u^1 (final velocity at time t)
    for i in range(d):
        for bc in bcs_u_corr:
            bc.apply(rhs_u_corr[i])
    for i in range(d):
        solve(lhs_u_corr[i], rhs_u_corr[i], u[1][i].vector()) # TODO: Setup solvers outside of loop

    print "done solving for u(t), t =", t
print "done"
