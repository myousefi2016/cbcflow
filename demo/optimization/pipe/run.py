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

# Configure scheme
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


# Configure problem
from pipe import Pipe

ppd = ParamDict(
    #dt = 1e-3,
    #T  = 1e-3 * 100,
    num_periods=0.002, #3.0,
    num_timesteps=2,#30,#100,
    alpha=1e-4,
    )
problem = Pipe(ppd)


# Configure postprocessing
pppd = ParamDict(casedir="pipe_results")
postprocessor = NSPostProcessor(pppd)

ppfield_pd = ParamDict(
    saveparams=ParamDict(
        save=False,#True,
        ),
    timeparams=ParamDict(
        step_frequency=1,#5,
        )
    )
#wss = WSS(params=ppfield_pd)
#velocity = Velocity(params=ppfield_pd)
#pressure = Pressure(params=ppfield_pd)

#postprocessor.add_fields([wss, velocity, pressure])
#postprocessor.add_fields([velocity, pressure])

if 0:
    parameters["optimization"]["test_gradient"] = True
    parameters["optimization"]["test_gradient_seed"] = 0.01

# Configure solver
npd = ParamDict(
    plot_solution=False,
    check_mem_usage=True,
    enable_annotation=enable_annotation,
    )
# Solve
solver = NSSolver(problem, scheme, postprocessor, params=npd)
sns = solver.solve()


# Try to replay
if 0 and enable_annotation:
    res = replay_dolfin()
    print "replay result =", res


# Optimization
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

if 1 and enable_annotation:
    # Fetch some stuff from scheme namespace
    V = sns["V"]
    Q = sns["Q"]
    controls = sns["controls"]
    state = sns["state"]
    t = sns["t"]

    u0, p0, p_out_coeffs = controls
    u, p = state
    d = len(u)

    J = Functional(problem.J(V, Q, t, u, p, controls))
    m = [InitialConditionParameter(u0c) for u0c in u0]
    m += [InitialConditionParameter(p0)]
    m += [ScalarParameter(p_coeff) for p_coeff in p_out_coeffs]

    # Dump annotation data
    if 1:
        adj_html("forward.html", "forward")
        adj_html("adjoint.html", "adjoint")

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
        Jm = Jred([u0c for u0c in u0] + [p0] + [p_coeff for p_coeff in p_out_coeffs])
        seed = 0.001
        dJdm = compute_gradient(J, m, forget=False)
        res = dolfin_adjoint.taylor_test(Jred, m, Jm, dJdm, seed=seed)
        print "taylor test result =", res

        dJdu0 = dJdm[:d]
        dJdp0 = dJdm[d]
        dJdp = dJdm[d+1:]
        print 'dJdu0 ='
        print dJdu0
        print 'dJdp0 ='
        print dJdp0
        print 'dJdp ='
        print dJdp
        print
        print 'dJdu0 ='
        print [norm(dj) for dj in dJdu0]
        print 'dJdp0 ='
        print [norm(dj) for dj in dJdp0]
        print 'dJdp ='
        print [norm(dj) for dj in dJdp]

    # Try to optimize
    if 1:
        Jred = ReducedFunctional(J, m,
                                 eval_cb=on_J_eval,
                                 derivative_cb=on_J_derivative,
                                 replay_cb=on_replay)
        m_opt = minimize(Jred, options={"disp":True}) # bounds=bounds, tol=1e-6,
        # TODO: Store m_opt through postprocessing framework instead of this adhoc implementation
        u0 = m_opt[:d]
        p0 = m_opt[d]
        p_coeffs = m_opt[d+1:]
        u0 = project(as_vector(u0), MixedFunctionSpace([V]*len(u0)))
        p = numpy.asarray([(t, sum(c*N for c,N in zip(p_out_coeffs, problem.pressure_basis(t))))
                            for t in numpy.linspace(problem.params.T0, problem.params.T, 100)])
        if not os.path.exists("optresults"):
            os.mkdir("optresults")
        File("optresults/u0.pvd") << u0
        File("optresults/p0.pvd") << p0
        open("optresults/pcoeffs.txt", "w").writelines(map(lambda x: "%g\n"%x, p_out_coeffs))
        numpy.savetxt("optresults/p.txt", p)
        p2 = numpy.loadtxt("optresults/p.txt")
        # TODO: Rerun forward problem with postprocessing to inspect final transient state

from headflow.core.utils import get_memory_usage
print "Memory usage at end of program:", get_memory_usage()

