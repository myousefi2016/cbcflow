#!/usr/bin/env python

from headflow import ParamDict, NSSolver, NSPostProcessor, all_schemes, show_problem
from headflow.postprocessing import *
import sys

all_scheme_names = sorted([s.__name__ for s in all_schemes])

def find_scheme(schemename):
    for s in all_schemes:
        if s.__name__ == schemename:
            return s
    return None

def solve_ns_problem(Problem, Scheme, args):
    # Input
    params = ParamDict(
        problem=Problem.default_params(),
        postproc=NSPostProcessor.default_params(),
        scheme=Scheme.default_params(),
        solver=NSSolver.default_params(),
        )
    #params.parse_cmdline(args) # FIXME: Parse commandline arguments

    # Setup
    problem = Problem(params.problem)
    scheme = Scheme(params.scheme)
    params.postproc.casedir = "results_demo_%s_%s" % (problem.shortname(), scheme.shortname())
    postproc = NSPostProcessor(params.postproc)
    solver = NSSolver(problem, scheme, postproc, params.solver)

    # FIXME: Need better way to configure what to show/save/compute for demos!

    # Add postprocessing fields (a few configurations)
    if 1:
        fp = dict(plot=True, save=True)#, start_timestep=1)
        #f = L2norm(VelocityError())
        #f = L2norm(Pressure(), fp)
        f = "t"
        postproc.add_fields([
            #f,
            #TimeIntegral(f, fp),
            TimeDerivative(f, fp),
            #SecondTimeDerivative(f, fp),
            ])
        print 'plan', postproc._plan
        print 'cache', postproc._cache
        print 'deps', postproc._dependencies
        print 'keys', postproc._sorted_fields_keys
        print 'fields', postproc._fields

    if 0:
        fp = dict(plot=True, save=True)
        postproc.add_fields([
            PressureError(),
            VelocityError(),
            L2norm("PressureError"),
            L2norm("VelocityError", fp),
            H1norm("VelocityError"),
            RunningMax("L2norm_PressureError", fp),
            RunningMax("H1norm_VelocityError", fp),
            DomainAvg("Pressure", fp),
            ])
    
    if 0:
        postproc.add_fields([
            Pressure(dict(plot=False, save=True)),
            Velocity(dict(plot=False, save=True)),
            ])

    if 0:
        postproc.add_fields([
            Pressure(dict(plot=True, save=True)),
            Velocity(dict(plot=True, save=True)),
            ])

    if 0:
        postproc.add_fields([
            Velocity(dict(save=True)),
            Pressure(dict(save=True)),
            H1norm("Velocity", dict(save=True)),
            Stress(dict(save=True)),
            WSS(dict(save=True)),
            Q(dict(save=True)),
            ])

    if 0:
        postproc.add_fields([
            Pressure(dict(plot=True)),
            Velocity(dict(plot=True)),
            VelocityCurl(dict(plot=True)),
            Lambda2(dict(plot=True)),
            ])

    if 0:
        postproc.add_fields([
            Pressure(dict(plot=True)),
            Velocity(dict(plot=True)),
            VelocityCurl(dict(plot=True)),
            VelocityDivergence(dict(plot=True)),
            L2norm("Velocity", dict(plot=True)),
            H1norm("Velocity", dict(plot=True)),
            H1seminorm("Velocity", dict(plot=True)),
            Strain(dict(plot=True)),
            Stress(dict(plot=True)),
            Q(dict(plot=True)),
            Delta(dict(plot=True)),
            Lambda2(dict(plot=True)),
            ])

    # Execution
    namespace = solver.solve()
    return namespace

def demo_main(Problem):
    args = sys.argv[1:]
    if not args:
        print "Usage:"
        print "  %s  <schemename>|show  <parameters>" % sys.argv[0]
        print "Valid schemes are:"
        print '\n'.join(all_scheme_names)
        return 1
    elif args[0] == "show":
        show_problem(Problem())
    else:
        Scheme = find_scheme(args[0])
        if Scheme is None:
            print "Couldn't find scheme %s." % args[0]
            return 2
        solve_ns_problem(Problem, Scheme, args[1:])
    return 0

if __name__ == "__main__":
    print "This file is not intended for running."
    # Usage: Simply call demo_main(Problem) in the demo problem file
