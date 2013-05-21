#!/usr/bin/env python

# Hack to run without installing, useful while working
import sys; sys.path.insert(0,"../../../site-packages")

import os
import numpy

enable_annotation = True

import dolfin
if enable_annotation:
    import dolfin_adjoint

from headflow import *
from headflow.dol import *
set_log_level(WARNING)
#set_log_level(PROGRESS)


from headflow.core.utils import get_memory_usage
print "Memory usage at top of program:", get_memory_usage()

parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True


# ====== Configure problem
from pipe import Problem
jp = ParamDict(
    alpha=1e-3,
    alpha_u_prior=1,
    alpha_u_grad=0,
    alpha_u_div=1e3,
    alpha_u_curl=0,
    alpha_u_wall=1,
    cyclic=0,
    )
ppd = ParamDict(
    #dt = 1e-3,
    #period=0.1,
    #num_periods=0.002, #3.0,
    num_timesteps=10,
    J=jp,
    pdim=3,
    )
problem = Problem(ppd)


# ====== Configure postprocessing
pppd = ParamDict(casedir="results_initial_run")
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


# ====== Configure postprocessing for final state
pppd = ParamDict(casedir="results_final_state")
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
postprocessor2.add_fields([velocity, pressure])


# ====== Configure scheme
spd = ParamDict()
scheme = CoupledPicard(spd)


if 0:
    parameters["optimization"]["test_gradient"] = True
    parameters["optimization"]["test_gradient_seed"] = 0.01


# ====== Configure solver
npd = ParamDict(
    plot_solution=False,
    check_mem_usage=True,
    enable_annotation=enable_annotation,
    )
# Solve
solver = NSSolver(problem, scheme, postprocessor, params=npd)
sns = solver.solve()

print "Memory usage after solve:", get_memory_usage()


# ====== Optimization

_J_eval_count = 0
_J_derivative_count = 0

def on_replay(var, func, m):
    print "/// Replay#%d at %s of %s, ||%s|| = %g" % (_J_eval_count, var.timestep, str(var), func.name(), norm(func))
    # TODO: Store func for timestep/iter

def on_J_eval(j, m):
    global _J_eval_count
    print "/// J evaluated #%d: %g" % (_J_eval_count, j)
    # TODO: Store m, j for iter
    _J_eval_count += 1

def on_J_derivative(j, dj, m):
    global _J_derivative_count
    def norm2(x):
        if isinstance(x, (float,int)):
            return abs(x)
        else:
            return norm(x)
    norms = map(lambda x: "%g" % norm2(x), dj)
    print "/// DJ evaluated #%d: %s" % (_J_derivative_count, norms)
    # TODO: Store m, j, dj for iter
    _J_derivative_count += 1

if enable_annotation:
    # Dump annotation data
    adj_html("forward.html", "forward")
    adj_html("adjoint.html", "adjoint")

if 0 and enable_annotation:
    # Try to replay
    res = replay_dolfin()
    print "replay result =", res

if 1 and enable_annotation:
    # Fetch some stuff from scheme namespace
    t = sns["t"]

    spaces = sns["spaces"]
    d = spaces.V.cell().d

    observations = sns["observations"]

    controls = sns["controls"]
    u0, p_out_coeffs = controls

    states = sns["states"]
    u, p = states

    J = Functional(problem.J(spaces, t, u, p, controls, observations))
    m = [InitialConditionParameter(u0c) for u0c in u0]
    m += [ScalarParameter(p_coeff) for p_coeff in p_out_coeffs]

    # Try to compute gradient
    if 0:
        #m = InitialConditionParameter(u0[0])
        m = ScalarParameter(p_out_coeffs[0])
        Jred = ReducedFunctional(J, m)
        dJdm = compute_gradient(J, m, forget=False)
        #Jm = Jred(m)
        #res = dolfin_adjoint.taylor_test(Jred, m, Jm, dJdm)#, seed=0.001)
        #print "taylor test result =", res

    if 0:
        Jred = ReducedFunctional(J, m)
        Jm = Jred([u0c for u0c in u0] + [p_coeff for p_coeff in p_out_coeffs])
        seed = 0.001
        dJdm = compute_gradient(J, m, forget=False)
        res = dolfin_adjoint.taylor_test(Jred, m, Jm, dJdm, seed=seed)
        print "taylor test result =", res

        dJdu0 = dJdm[:d]
        dJdp = dJdm[d:]
        print 'dJdu0 ='
        print dJdu0
        print 'dJdp ='
        print dJdp
        print
        print 'dJdu0 ='
        print [norm(dj) for dj in dJdu0]
        print 'dJdp ='
        print [norm(dj) for dj in dJdp]

    # Try to replay through functional
    if 0:
        Jred = ReducedFunctional(J, m,
                                 eval_cb=on_J_eval,
                                 derivative_cb=on_J_derivative,
                                 replay_cb=on_replay)
        Jm = Jred([u0c for u0c in u0] + [p_coeff for p_coeff in p_out_coeffs])
        print "Jm =", Jm
        crash

    # Try to optimize
    if 1:
        Jred = ReducedFunctional(J, m,
                                 eval_cb=on_J_eval,
                                 derivative_cb=on_J_derivative,
                                 replay_cb=on_replay)
        m_opt = minimize(Jred, options={"disp":True}) # bounds=bounds, tol=1e-6,

    if 1:
        # TODO: Store m_opt through postprocessing framework instead of this adhoc implementation
        u0 = m_opt[:d]
        p_coeffs = m_opt[d:]
        u0 = project(as_vector(u0), spaces.V)
        p = numpy.asarray([(t, sum(c*N for c,N in zip(p_out_coeffs, problem.pressure_basis(t))))
                            for t in numpy.linspace(problem.params.T0, problem.params.T, 100)])
        if not os.path.exists("results_control"):
            os.mkdir("results_control")
        File("results_control/u0.pvd") << u0
        open("results_control/pcoeffs.txt", "w").writelines(map(lambda x: "%g\n"%x, p_out_coeffs))
        numpy.savetxt("results_control/p.txt", p)
        p2 = numpy.loadtxt("results_control/p.txt")

    if 1:
        # TODO: Rerun forward problem with postprocessing to inspect final transient state

        problem.set_controls(m_opt)

        # ====== Configure solver
        npd = ParamDict(
            plot_solution=False,
            check_mem_usage=True,
            enable_annotation=False,
            )
        # Solve
        solver = NSSolver(problem, scheme, postprocessor2, params=npd)
        sns = solver.solve()

        print "Memory usage after solve:", get_memory_usage()

print "Memory usage at end of program:", get_memory_usage()
