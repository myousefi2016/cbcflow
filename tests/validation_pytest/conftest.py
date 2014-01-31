from cbcflow import *
import sys
import pytest
import itertools
from collections import defaultdict
sys.path.insert(0, "../../demo/undocumented/Poiseuille2D")
sys.path.insert(0, "../../demo/undocumented/Poiseuille3D")
sys.path.insert(0, "../../demo/undocumented/Womersley2D")
sys.path.insert(0, "../../demo/undocumented/Womersley3D")
sys.path.insert(0, "../../demo/undocumented/Beltrami")
sys.path.insert(0, "../../demo/undocumented/LidDrivenCavity")
from poiseuille2d import Poiseuille2D
from poiseuille3d import Poiseuille3D
from womersley2d import Womersley2D
from womersley3d import Womersley3D
from beltrami import Beltrami
from liddrivencavity import LidDrivenCavity


# TODO: Move to utils, this could be useful elsewhere
def parse_parameterized(option, opt_str, value, parser):
    parameterized = {}
    nargs = 0
    # Start parsing arguments
    for arg in parser.rargs:
        # Stop parsing if we've reached a new option
        if arg[0] == '-':
            break
        nargs += 1
        
        # Create new key if this is not recognized as a parameter
        if "=" not in arg:
            key = arg
            parameterized[key] = ParamDict()
            continue
        
        # Add parameters
        param, value = arg.split('=')
        if value.isdigit():
            parameterized[key][param] = int(value)
        else:
            try:
                parameterized[key][param] = float(value)
            except ValueError:
                parameterized[key][param] = value
    
    # Removed parsed args from parser
    del parser.rargs[:nargs]
    setattr(parser.values, option.dest, parameterized)

def pytest_addoption(parser):
    parser.addoption("--type", action="store", default="changed", 
        choices=["changed", "update"], dest="type", help="check if reference has changed, or update reference")
    parser.addoption("--testsize", action="store", default="all",
        choices=["all", "fast", "slow", "slowest", "debug"], dest="testsize", help="run all, slow or fast test")
    
    # TODO: Add option for validation mode
    #parser.addoption("--mode", action="store", default="regression", choices=["regression", "validation"],
    #                 dest="mode", help="Run regression test or generate validation figures")
    
    parser.addoption("--schemes", action="callback", callback=parse_parameterized, dest="schemes")
    parser.addoption("--problems", action="callback", callback=parse_parameterized, dest="problems")
    
def updated_paramdict(params, add_dict):
    params.update(add_dict)
    return params

def make_factory(cl, **kwargs):
    def factory():
        return cl(**kwargs)
    return factory


def create_problem_factories(cmdline_problem_args):   
    if not cmdline_problem_args:
        # Default problem factories
        problem_factories = [
            lambda refinement_level,dt: LidDrivenCavity(ParamDict(refinement_level=refinement_level, dt=dt, T=2.5, mu=1./1000., rho=1.0)),
            lambda refinement_level,dt: Poiseuille2D(ParamDict(refinement_level=refinement_level, dt=dt, T=None, num_periods=0.2)),
            lambda refinement_level,dt: Poiseuille3D(ParamDict(refinement_level=refinement_level, dt=dt, T=None, num_periods=0.2)),
            lambda refinement_level,dt: Womersley2D(ParamDict(refinement_level=refinement_level, dt=dt, T=None, num_periods=1.0)),
            lambda refinement_level,dt: Womersley3D(ParamDict(refinement_level=refinement_level, dt=dt, T=None, num_periods=1.0)), 
            lambda refinement_level,dt: Beltrami(ParamDict(refinement_level=refinement_level, dt=dt, T=1.0)),
        ]
    else:
        # Create factories from what has been parsed from --problems option at command line
        problem_factories = []
        for problem_name, params in cmdline_problem_args.items():
            Problem = eval(problem_name)
            factory = (lambda Problem, params: lambda refinement_level, dt: Problem(updated_paramdict(params, dict(refinement_level=refinement_level, dt=dt))))(Problem, params)
            problem_factories.append(factory)

        
    return problem_factories


def create_scheme_factories(cmdline_scheme_args):
    if not cmdline_scheme_args:
        # Default schemes run
        scheme_factories = [
            lambda: IPCS(),
            lambda: Yosida(),
            lambda: IPCS_Stable(ParamDict(theta=1.0)),
            lambda: IPCS_Stable(ParamDict(theta=0.5)),
        ]
    else:
        scheme_factories = []

        for scheme_name, params in cmdline_scheme_args.items():
            Scheme = eval(scheme_name)
            #print params
            factory = (lambda Scheme, params : lambda: Scheme(params))(Scheme, params)
            scheme_factories.append(factory)
        
    return scheme_factories



def pytest_generate_tests(metafunc):
    if metafunc.function.func_name != "test_run_problem":
        return
    else:
        test_run_problem = metafunc

    # Create problem and scheme factoriehe pos
    problem_factories = create_problem_factories(metafunc.config.option.problems)
    scheme_factories = create_scheme_factories(metafunc.config.option.schemes)

    testsize = metafunc.config.option.testsize    
    if testsize == "all":
        params = dict(
            refinements = range(5),
            dts = [0.1, 0.05, 0.025, 0.0125, 0.00625],
        )
    elif testsize == "fast":
        params = dict(
            refinements = range(3),
            dts = [0.1, 0.05],
        )
    elif testsize == "slow":
        params = dict(
            refinements = range(3,5),
            dts = [0.025, 0.0125, 0.00625],
        )
    elif testsize == "debug":
        params = dict(
            refinements = [0],
            dts = [0.1]
        )


    test_run_problem.parametrize("scheme_factory", scheme_factories)
    test_run_problem.parametrize("problem_factory", problem_factories)
    test_run_problem.parametrize("refinement_level", params["refinements"])
    test_run_problem.parametrize("dt", params["dts"])
    
    
    