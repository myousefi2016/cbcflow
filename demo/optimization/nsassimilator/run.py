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
        # Optimization process params
        minimize_tolerance=1e-6,

        # Pick alternative actions, skipping optimization:
        replay_only=False,
        evaluate_only=False,
        enable_test_gradient=False,

        verbose=False,
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

    # Interpret controls
    d = spaces.d
    u0 = controls[:d]
    p_coeffs = controls[d:]

    # Store velocity
    for k in spaces.dims:
        File(os.path.join(casedir, "u0_%d.xml.gz" % k)) << u0[k]
    u0vec = project(as_vector(u0), spaces.V)
    f = File(os.path.join(casedir, "u0.pvd"))
    f << u0vec

    # Store pressure coeffs
    f = open(os.path.join(casedir, "p_coeffs.txt", "w"))
    f.writelines(map(lambda x: "%g\n"%x, p_coeffs))

    # Store pressure evaluated at timesteps
    p = numpy.asarray([(t, sum(c*N for c,N in zip(p_coeffs, problem.pressure_basis(t))))
                        for t in timesteps])
    numpy.savetxt(os.path.join(casedir, "p.txt", p))

    #p2 = numpy.loadtxt("results_control/p.txt")


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
    u0, p_out_coeffs = controls

    states = sns["states"]
    u, p = states


    # ====== Optimization callbacks
    _J_eval_count = 0
    _J_derivative_count = 0

    def on_replay(var, func, m):
        fn = norm2(func)
        headflow_print("/// Replay#%d at %s of %s, ||%s|| = %g" % (_J_eval_count, var.timestep, str(var), func.name(), fn))
        # TODO: Use postprocessing framework to store subset of functions for subset of timestep/iter

    def on_J_eval(j, m):
        headflow_print("/// J evaluated #%d: %g" % (_J_eval_count, j))
        casedir = os.path.join(controls_output_dir, "iteration%d" % _J_eval_count)
        store_controls(controls, j, spaces, timesteps, problem, casedir)
        _J_eval_count += 1

    def on_J_derivative(j, dj, m):
        # FIXME: Store dj!
        norms = map(lambda x: norm2(x), dj)
        norms = map(lambda x: "%g" % x, norms)
        headflow_print("/// DJ evaluated #%d: %s" % (_J_derivative_count, norms))
        _J_derivative_count += 1


    # ====== Setup reduced functional
    J = Functional(problem.J(spaces, t, u, p, controls, observations))
    m = [InitialConditionParameter(u0c) for u0c in u0]
    m += [ScalarParameter(p_coeff) for p_coeff in p_out_coeffs]
    Jred = ReducedFunctional(J, m,
                             eval_cb=on_J_eval,
                             derivative_cb=on_J_derivative,
                             replay_cb=on_replay)


    # ====== Test replay through evaluation of reduced functional
    if params.assimilator.evaluate_only:
        Jm = Jred([u0c for u0c in u0] + [p_coeff for p_coeff in p_out_coeffs])
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
    p.scale=1e4,
    p.pdim=1,
    sys.exit(run(p))
