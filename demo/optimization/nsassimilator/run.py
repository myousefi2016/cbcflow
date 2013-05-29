#!/usr/bin/env python

import os, sys, numpy, scipy


# ====== Import dolfin, dolfin_adjoint and headflow in proper order (order is important!)

# Hack to run without installing, useful while working
#import sys; sys.path.insert(0,"../../../site-packages")

import dolfin
import dolfin_adjoint
from headflow import *
from headflow.dol import *
from headflow.core.utils import get_memory_usage, headflow_print


# ====== Import problem and scheme
from problem import Problem # This limits us to a single problem spec in a single directory
from headflow import CoupledPicard as Scheme


# ====== Collection of default parameters
def default_params():

    # Configure solver
    npd = NSSolver.default_params()
    npd.plot_solution = False
    npd.check_mem_usage = True
    npd.enable_annotation = True

    # Configure top level behaviour
    pa = ParamDict(
        casedir="default_casedir",
        verbose=False,

        # Optimization process params
        minimize_tolerance=1e-6,

        # Pick alternative actions, skipping optimization:
        replay_only=False,
        evaluate_only=False,
        enable_test_gradient=False,
        )

    # Collect params
    p = ParamDict(
        problem=Problem.default_params(),
        scheme=Scheme.default_params(),
        solver=npd,
        assimilator=pa,
        )

    return p

def ensure_dir(casedir):
    if os.path.exists(casedir):
        print "WARNING! Directory %s already exists." % casedir
    else:
        os.mkdir(casedir)

def norm2(x):
    if isinstance(x, (float,int)):
        return abs(x)
    else:
        return norm(x)

# TODO: Store m_opt through postprocessing framework instead of this adhoc implementation
def store_controls(controls, spaces, timesteps, problem, casedir):
    ensure_dir(casedir)

    def fn(name):
        return os.path.join(casedir, name)

    def save_array(name, data):
        numpy.savetxt(fn(name), data)

    def _assemble(form):
        return assemble(form, mesh=problem.mesh, annotate=False)

    def _project(expr, space):
        return project(expr, space, annotate=False)

    n = FacetNormal(problem.mesh)
    ds = problem.ds

    # Get some dimensions
    d = spaces.d
    num_p_controls = len(problem.control_boundaries)
    pdim = problem.params.pdim

    # Interpret controls (flat array)
    assert len(controls) == d + pdim*num_p_controls
    u0 = as_vector(controls[:d])
    p_coeffs = [as_vector(controls[d+k*pdim: d+(k+1)*pdim])
                for k in xrange(num_p_controls)]

    # Store velocity components # TODO: Use xdmf?
    for k in spaces.dims:
        File(fn("u0_%d.xml.gz" % k)) << u0[k]

    # Store velocity vector for visualization
    File(fn("u0.pvd")) << _project(u0, spaces.V_CG1)

    # Store some postprocessed quantities for convenient inspection
    File(fn("div_u0.pvd")) << _project(div(u0), spaces.U_CG1)

    # Store norms
    with open(fn("norms.txt"), "w") as f:
        u0norm = sqrt(_assemble(u0**2*dx()))
        f.write("||u0||     = %g\n" % u0norm)
        u0divnorm = sqrt(_assemble(div(u0)**2*dx()))
        f.write("||div u0|| = %g\n" % u0divnorm)
        u0n = _assemble(dot(u0,n)*ds())
        f.write("int_Gamma (u0 . n) = %g\n" % u0n)

    # Store flow rates
    with open(fn("flowrates.txt"), "w") as f:
        for name, boundaries in [("walls", problem.wall_boundaries),
                                 ("given", problem.given_pressure_boundaries),
                                 ("controls", problem.control_boundaries)]:
            f.write("%s\n" % name)
            for r in boundaries:
                Qr = _assemble(dot(u0, n)*ds(r))
                f.write("%d  %.6e\n" % (r, Qr))

    # Store raw pressure coeffs
    pc = numpy.zeros((num_p_controls, pdim))
    for k in xrange(num_p_controls):
        for i in xrange(pdim):
            pc[k,i] = p_coeffs[k][i]
    save_array("p_coeffs.txt", pc)

    # Store pressure and basis functions evaluated at timesteps
    pv = numpy.zeros((len(timesteps), 1+num_p_controls))
    Nv = numpy.zeros((len(timesteps), 1+pdim))
    for j, t in enumerate(timesteps):
        Ns = problem.pressure_basis(t)
        pv[j,0] = t
        Nv[j,0] = t
        for i in xrange(pdim):
            Nv[j,i+1] = Ns[i]
        for k in xrange(num_p_controls):
            pk = sum(c*N for c,N in zip(p_coeffs[k],Ns))
            pv[j,k+1] = pk
    save_array("p_controls.txt", pv)
    save_array("p_basis.txt", Nv)

    # Debugging plots TODO: Move to analysis script for inspecting
    if 0:
        pv = numpy.loadtxt(fn("p_controls.txt"))
        Nv = numpy.loadtxt(fn("p_basis.txt"))
        assert pv.shape == (len(timesteps), 1+num_p_controls)
        assert Nv.shape == (len(timesteps), 1+pdim)
        import pylab
        pylab.figure(1)
        pylab.hold(True)
        for k in xrange(num_p_controls):
            pylab.plot(pv[0,:], pv[1+k,:])
        pylab.legend(["p %d" % k for k in xrange(num_p_controls)])
        pylab.hold(False)
        pylab.figure(2)
        pylab.hold(True)
        for k in xrange(pdim):
            pylab.plot(Nv[0,:], Nv[1+k,:])
        pylab.legend(["N %d" % k for k in xrange(pdim)])
        pylab.hold(False)

# ====== The overall driver function
def run(params):

    # ====== Set some global parameters
    if params.assimilator.verbose:
        set_log_level(PROGRESS)
    else:
        set_log_level(WARNING)

    parameters["form_compiler"]["optimize"]     = True
    parameters["form_compiler"]["cpp_optimize"] = True
    parameters["form_compiler"]["cpp_optimize_flags"] = "-O3 -march=native -fno-math-errno"


    # ====== Enable (costly) taylor test of gradient in optimization algorithm
    if params.assimilator.enable_test_gradient:
        parameters["optimization"]["test_gradient"] = True
        parameters["optimization"]["test_gradient_seed"] = 0.001


    # ====== Setup top level casedir
    casedir = params.assimilator.casedir
    controls_output_dir = os.path.join(casedir, "controls")
    ensure_dir(casedir)
    ensure_dir(controls_output_dir)


    # ====== Configure postprocessing for initial forward run
    pppd = ParamDict(
        casedir = os.path.join(casedir, "initial"),
        )
    postprocessor = NSPostProcessor(pppd)
    ppfield_pd = ParamDict(
        saveparams=ParamDict(
            save=True,
            ),
        timeparams=ParamDict(
            step_frequency=1,
            )
        )
    velocity = Velocity(params=ppfield_pd)
    pressure = Pressure(params=ppfield_pd)
    postprocessor.add_fields([velocity, pressure])


    # ====== Solve forward problem with annotation
    problem = Problem(params.problem)
    scheme = Scheme(params.scheme)
    params.solver.enable_annotation = True
    solver = NSSolver(problem, scheme, postprocessor, params=params.solver)
    sns = solver.solve()

    # Always dump annotation data for debugging
    adj_html("forward.html", "forward")
    adj_html("adjoint.html", "adjoint")


    # ====== Test replay functionality
    if params.assimilator.replay_only:
        res = replay_dolfin()
        headflow_print("replay result =", res)
        return 0


    # ====== Fetch variables from scheme namespace (slightly hacky solution)
    t = sns["t"]

    spaces = sns["spaces"]
    d = spaces.V.cell().d

    timesteps = sns["timesteps"]

    observations = sns["observations"]

    controls = sns["controls"]
    u0, p_coeffs = controls
    num_p_controls = len(problem.control_boundaries)
    pdim = problem.params.pdim

    states = sns["states"]
    u, p = states


    # ====== Optimization callbacks
    def on_replay(var, func, m):
        fn = norm2(func)
        headflow_print("/// Replay#%d at %s of %s, ||%s|| = %g" % (on_J_eval.num_calls, var.timestep, str(var), func.name(), fn))
        # TODO: Use postprocessing framework to store subset of functions for subset of timestep/iter

    def on_J_eval(j, m):
        # FIXME: Store j?
        c = on_J_eval.num_calls
        headflow_print("/// J evaluated #%d: %g" % (c, j))
        casedir = os.path.join(controls_output_dir, "J_call_%d" % c)
        store_controls(m, spaces, timesteps, problem, casedir)
        on_J_eval.num_calls += 1
    on_J_eval.num_calls = 0

    def on_J_derivative(j, dj, m):
        # FIXME: Store dj?
        norms = map(lambda x: norm2(x), dj)
        norms = map(lambda x: "%g" % x, norms)
        headflow_print("/// DJ evaluated #%d: %s" % (on_J_derivative.num_calls, norms))
        on_J_derivative.num_calls += 1
    on_J_derivative.num_calls = 0


    # ====== Setup reduced functional
    J = Functional(problem.J(spaces, t, u, p, controls, observations))
    m = [InitialConditionParameter(u0c) for u0c in u0]
    for k in xrange(num_p_controls):
        m += [ScalarParameter(p_coeffs[k][i]) for i in xrange(pdim)]

    Jred = ReducedFunctional(J, m,
                             eval_cb=on_J_eval,
                             derivative_cb=on_J_derivative,
                             replay_cb=on_replay)


    # ====== Test replay through evaluation of reduced functional
    if params.assimilator.evaluate_only:
        mval = list(u0)
        for k in xrange(num_p_controls):
            mval += list(p_coeffs[k])
        Jm = Jred(mval)
        headflow_print("Jm = %g" % Jm)
        return 0


    # ====== Execute optimization problem
    m_opt = minimize(Jred,
                     tol=params.assimilator.minimize_tolerance,
                     options={"disp":True},
                     ) # bounds=bounds,


    # ====== Store final control variables
    final_casedir = os.path.join(controls_output_dir, "final")
    store_controls(m_opt, spaces, timesteps, problem, final_casedir)


    # ====== Configure postprocessing for final state
    pppd = ParamDict(
        casedir = os.path.join(casedir, "final"),
        )
    postprocessor2 = NSPostProcessor(pppd)
    ppfield_pd = ParamDict(
        saveparams=ParamDict(
            save=True,
            ),
        timeparams=ParamDict(
            step_frequency=1,
            )
        )
    velocity = Velocity(params=ppfield_pd)
    pressure = Pressure(params=ppfield_pd)
    postprocessor.add_fields([velocity, pressure])


    # ====== Rerun forward problem with postprocessing to inspect final transient state
    problem.set_controls(m_opt)
    params.solver.enable_annotation = False
    solver = NSSolver(problem, scheme, postprocessor, params=params.solver)
    sns = solver.solve()


    # ====== End of program
    print "Memory usage at end of program:", get_memory_usage()
    return 0

if __name__ == "__main__":

    p = default_params()

    p.problem.num_timesteps = 1
    p.problem.scale = 1e4
    p.problem.pdim = 1

    p.assimilator.casedir = "single_timestep_test"
    p.assimilator.minimize_tolerance = 1e-4

    sys.exit(run(p))
