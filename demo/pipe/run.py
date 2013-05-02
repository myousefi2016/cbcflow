#!/usr/bin/env python

# Hack to run without installing, useful while working
import sys; sys.path.insert(0,"../../site-packages")

enable_annotation = True

import dolfin
if enable_annotation:
    import dolfin_adjoint

from headflow import *
from headflow.dol import *
set_log_level(100)


# Configure scheme
spd = ParamDict(
    u_degree=1,
    solver_p=("lu", "default"),
    #solver_p=("gmres", "ilu"),
    #solver_p=("cg", "amg"),
    #solver_p=("gmres", "amg"),
    )
scheme = PenaltyIPCS(spd)


# Configure problem
from pipe import Pipe

ppd = ParamDict(
    #dt = 1e-3,
    #T  = 1e-3 * 100,
    num_periods=0.01,#3.0,
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
    rep = replay_dolfin()
    print "rep =", rep


# Optimization
def update_model_eval(*args):
    print "update_model_eval: ", args

def update_derivative_eval(*args):
    print "update_derivative_eval: ", args

if 1 and enable_annotation:
    # Fetch some stuff from scheme namespace
    V = sns["V"]
    Q = sns["Q"]
    controls = sns["controls"]
    state = sns["state"]
    t = sns["t"]

    u0, p_out_coeffs = controls
    u, p = state

    J = Functional(problem.J(V, Q, t, u, p, controls))
    m = [InitialConditionParameter(u0c) for u0c in u0]
    m += [InitialConditionParameter(p_coeff) for p_coeff in p_out_coeffs]

    # Try to compute gradient
    if 1:
        dJdm = compute_gradient(J, m) # FIXME: This fails in dolfin-adjoint
        dJdu = dJdm[:d]
        dJdp = dJdm[d:]
        print 'dJdu ='
        print [norm(dj.vector()) for dj in dJdu]
        print 'dJdp ='
        print [norm(dj.vector()) for dj in dJdp]

    # Try to optimize
    if 0:
        Jred = ReducedFunctional(J, m,
                                 eval_cb=update_model_eval,
                                 eval_derivative_cb=update_derivative_eval)
        m_opt = minimize(Jred, options={"disp":True}) # bounds=bounds, tol=1e-6,
