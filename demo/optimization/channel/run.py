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
from channel import Problem
jp = ParamDict(
    alpha=1e-4,
    cyclic=1,
    )
ppd = ParamDict(
    dt = 1e-3,
    period=0.1,
    #num_periods=0.002, #3.0,
    num_timesteps=100,
    J=jp,
    )
problem = Problem(ppd)


# ====== Configure postprocessing for initial state
pppd = ParamDict(casedir="results_initial_state")
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
spd = ParamDict(
    u_degree=1,

    #solver_u_tent=("lu", "default"),
    #solver_u_corr=("lu", "default"),
    #solver_p=("lu", "default"),

    solver_u_tent=("gmres", "ilu"),
    solver_u_corr=("gmres", "ilu"),
    solver_p=("lu", "default"),

    #solver_u_tent=("gmres", "amg"),
    #solver_u_corr=("gmres", "amg"),
    #solver_p=("lu", "default"),

    #solver_p=("gmres", "default"),
    #solver_p=("cg", "default"),
    #solver_p=("gmres", "ilu"),
    #solver_p=("cg", "amg"),
    #solver_p=("gmres", "amg"),
    )
#scheme = PenaltyIPCS(spd)
scheme = SegregatedPenaltyIPCS(spd)


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

def on_replay(var, func, m):
    print "/// Replay %s: %s, %s (%g)" % (var.timestep, str(var), func.name(), norm(func))
    # TODO: Store func for timestep/iter

def on_J_eval(j, m):
    print "/// J evaluated: %g" % j
    # TODO: Store m, j for iter

def on_J_derivative(j, dj, m):
    def norm2(x):
        if isinstance(x, (float,int)):
            return abs(x)
        else:
            return norm(x)
    print "/// DJ evaluated: %s" % map(lambda x: "%g" % norm2(x), dj)
    # TODO: Store m, j, dj for iter

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
    V = sns["V"]
    Q = sns["Q"]
    controls = sns["controls"]
    state = sns["state"]
    t = sns["t"]

    u0, p_out_coeffs = controls
    u, p = state
    d = len(u)

    J = Functional(problem.J(V, Q, t, u, p, controls))
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
        u0 = project(as_vector(u0), MixedFunctionSpace([V]*d))
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

