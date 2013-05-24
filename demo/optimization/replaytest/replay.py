import ufl
from dolfin import *
from dolfin_adjoint import *

# This is a trick to handle automatic timestep annotation
def Time(t0=0.0):
    t = Constant(t0, name="TIME")
    t._prev_value = t0
    t._assigned_to = False
    return t

def assign_time(t, tvalue):
    if t._assigned_to:
        # Annotate previous timestep is done
        adj_inc_timestep(t._prev_value)
    else:
        # Annotate the beginning of time
        t._assigned_to = True
        adj_start_timestep(t._prev_value)
    # Update time constant to reflect modern times
    t.assign(tvalue)
    t._prev_value = tvalue

def finalize_time(t):
    # Make sure we have annotated the beginning of time
    if not t._assigned_to:
        t._assigned_to = True
        adj_start_timestep(t._prev_value)
    # Annotate the end-time is here
    adj_inc_timestep(t._prev_value, finished=True)
    # Time constant needs no updating anymore
    t._prev_value = None

def main():
    mesh = Mesh("pipe_0.2.xml.gz")
    facet_domains = mesh.domains().facet_domains()
    ds = ufl.ds[facet_domains]
    n  = FacetNormal(mesh)

    t = Time(t0=0.0)

    d = mesh.ufl_cell().geometric_dimension()
    dims = range(d)
    U = FunctionSpace(mesh, "CG", 2)
    V = VectorFunctionSpace(mesh, "CG", 2)
    Q = FunctionSpace(mesh, "CG", 1)
    W = V*Q

    # Test and trial functions
    u, p = TrialFunctions(W)
    v, q = TestFunctions(W)

    # Solution functions
    up0 = Function(W, name="up0") # Previous timestep
    up1 = Function(W, name="up1") # Last iterate of fixed-point in current timestep
    u0, p0 = split(up0)
    u1, p1 = split(up1)

    # Apply initial conditions and use it as initial guess
    uinit = [Function(U, name="ui_%d"%i) for i in xrange(d)]
    pinit = Function(Q, name="pinit")
    upinit = project(as_vector(uinit + [pinit]), W)
    up0.assign(upinit)
    up1.assign(up0)

    # Define boundary conditions
    c0 = Constant(0.0)
    p_out = Constant(0.0, name="p_out")
    p_out.assign(-5.6)
    bcu = [DirichletBC(W.sub(0).sub(d), c0, facet_domains, 0) for d in dims]
    Lbc = -dot(p_out*n, v)*ds(1)

    # Problem parameters
    nu = Constant(0.035, name="nu")
    k  = Constant(0.1, name="dt")

    # Picard linearization of Navier-Stokes, F = a*u - L = 0
    eqchoice = 1
    if eqchoice == 1:
        # Not scaled by k
        a = (
            dot((1.0/k)*u + (grad(u)*u1), v)*dx
            + nu*inner(grad(u), grad(v))*dx
            - p*div(v)*dx - q*div(u)*dx
            )
        L = dot((1.0/k)*u0, v)*dx + Lbc
    if eqchoice == 2:
        # Scaled by k (smaller residual, nonlinear solver hits absolute stopping criteria faster)
        a = (
            dot(u + k*(grad(u)*u1), v)*dx
            + (k*nu)*inner(grad(u), grad(v))*dx
            - (k*p)*div(v)*dx - (k*div(u))*q*dx
            )
        L = dot(u0, v)*dx + k*Lbc
    if eqchoice == 3:
        # Stokes
        a = (
            + nu*inner(grad(u), grad(v))*dx
            - p*div(v)*dx - q*div(u)*dx
            )
        L = Lbc
    F = action(a, up1) - L

    form_compiler_parameters = {
            "optimize": True,
            "cpp_optimize": True,
            "cpp_optimize_flags": "-O3 -march=native -fno-math-errno -ffinite-math-only -fno-rounding-math -fno-signaling-nans",
            "quadrature_degree": "auto",
            }

    picard_problem = NonlinearVariationalProblem(F, up1, bcu, J=a,
                                                 form_compiler_parameters=form_compiler_parameters)
    solver = NonlinearVariationalSolver(picard_problem)

    solver.parameters["lu_solver"]["report"] = True
    solver.parameters["lu_solver"]["verbose"] = True
    solver.parameters["newton_solver"]["report"] = True

    solver.parameters["lu_solver"]["reuse_factorization"] = True
    solver.parameters["lu_solver"]["same_nonzero_pattern"] = True
    solver.parameters["reset_jacobian"] = False

    solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-5
    solver.parameters["newton_solver"]["relative_tolerance"] = 1e-5

    for timestep in range(1,2):
        assign_time(t, timestep*float(k))
        solver.solve()
        up0.assign(up1)
    finalize_time(t)

    controls = uinit, p_out
    return u0, controls

# Run
u, controls = main()
adj_html("forward.html", "forward")
adj_html("adjoint.html", "adjoint")

# Replay
#res = replay_dolfin(forget=False)
#print "replay result =", res

# Test gradient
parameters["optimization"]["test_gradient"] = True
parameters["optimization"]["test_gradient_seed"] = 0.001 # TODO: Test smaller seed to see convergence results

uinit, p_out = controls
uinit = as_vector(uinit)
J = Functional(u**2*dx*dt + uinit**2*dx*dt[START_TIME] + p_out**2*ds(1)*dt[START_TIME])
m = [InitialConditionParameter(u0c) for u0c in uinit]
m += [ScalarParameter(p_coeff) for p_coeff in [p_out]]

Jred = ReducedFunctional(J, m)
m_opt = minimize(Jred, options={"disp":True})
