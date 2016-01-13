
import pytest

import sys
import os
import itertools
from collections import defaultdict
import argparse

from cbcflow import *
from cbcpost import ParamDict


# Hacks to import problem classes from demo directories
demo_dir = os.path.abspath(os.path.join(os.path.split(__file__)[0], "..", "..", "demo"))

# Documented demos:
ddemo_dir = os.path.join(demo_dir, "documented")
sys.path.insert(0, os.path.join(ddemo_dir, "Womersley3D"))

# Undocumented demos:
udemo_dir = os.path.join(demo_dir, "undocumented")
sys.path.insert(0, os.path.join(udemo_dir, "Poiseuille2D"))
sys.path.insert(0, os.path.join(udemo_dir, "Poiseuille3D"))
sys.path.insert(0, os.path.join(udemo_dir, "Womersley2D"))
sys.path.insert(0, os.path.join(udemo_dir, "Beltrami"))
sys.path.insert(0, os.path.join(udemo_dir, "LidDrivenCavity"))

from Poiseuille2D import Poiseuille2D
from Poiseuille3D import Poiseuille3D
from Womersley2D import Womersley2D
from Womersley3D import Womersley3D
from Beltrami import Beltrami
from LidDrivenCavity import LidDrivenCavity


# TODO: Move to utils, this could be useful elsewhere
#def parse_parameterized(option, opt_str, value, parser):
#def parse_parameterized(parser, dest=None, option_strings=None):

class ParseParameterized(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        parameterized = {}
        # Start parsing arguments
        for arg in values:
            # Stop parsing if we've reached a new option
            #if arg[0] == '-':
            #    break
            #nargs += 1

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
        setattr(namespace, self.dest, parameterized)


def pytest_addoption(parser):
    parser.addoption("--testsize", action="store", default="fast",
        choices=["all", "fast", "slow", "slowest", "debug"], dest="testsize", help="run all, slow or fast test")

    parser.addoption("--schemes", default=None, nargs="*", action=ParseParameterized, dest="schemes")
    parser.addoption("--problems", default=None, nargs="*", action=ParseParameterized, dest="problems")

def create_default_problem_factories():
    problem_factories = [
            #lambda refinement_level,dt: LidDrivenCavity(ParamDict(refinement_level=refinement_level, dt=dt, T=2.5, mu=1./1000., rho=1.0)),
            lambda refinement_level,dt: Poiseuille2D(ParamDict(refinement_level=refinement_level, dt=dt, T=None, num_periods=0.2)),
            lambda refinement_level,dt: Poiseuille3D(ParamDict(refinement_level=refinement_level, dt=dt, T=None, num_periods=0.2)),
            lambda refinement_level,dt: Womersley2D(ParamDict(refinement_level=refinement_level, dt=dt, T=None, num_periods=1.0)),
            lambda refinement_level,dt: Womersley3D(ParamDict(refinement_level=refinement_level, dt=dt, T=None, num_periods=1.0)),
            lambda refinement_level,dt: Beltrami(ParamDict(refinement_level=refinement_level, dt=dt, T=1.0)),
            ]
    return problem_factories

def create_custom_problem_factories(cmdline_problem_args):
    # Create factories from what has been parsed from --problems option at command line
    problem_factories = []
    for problem_name, params in cmdline_problem_args.items():
        Problem = eval(problem_name)
        factory = (lambda Problem, params:
                       (lambda refinement_level, dt:
                        Problem(params.update(refinement_level=refinement_level, dt=dt)))
                   )(Problem, params)
        problem_factories.append(factory)
    return problem_factories

def create_problem_factories(cmdline_problem_args):
    if cmdline_problem_args:
        return create_custom_problem_factories(cmdline_problem_args)
    else:
        return create_default_problem_factories()

def create_default_scheme_factories():
    scheme_factories = [
        lambda: IPCS_Naive(),
        #lambda: Yosida(),
        lambda: IPCS(ParamDict(theta=1.0)),
        lambda: IPCS(ParamDict(theta=0.5)),
        ]
    return scheme_factories

def create_custom_scheme_factories(cmdline_scheme_args):
    scheme_factories = []
    for scheme_name, params in cmdline_scheme_args.items():
        Scheme = eval(scheme_name)
        factory = (lambda Scheme, params:
                   (lambda: Scheme(params))
                   )(Scheme, params)
        scheme_factories.append(factory)
    return scheme_factories

def create_scheme_factories(cmdline_scheme_args):
    if cmdline_scheme_args:
        return create_custom_scheme_factories(cmdline_scheme_args)
    else:
        return create_default_scheme_factories()

def pytest_generate_tests(metafunc):
    if metafunc.function.func_name != "test_run_problem":
        return
    else:
        test_run_problem = metafunc

    # Create problem and scheme factories
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
    else:
        raise ValueError("Does not recognize testsize=%s." %testsize)

    test_run_problem.parametrize("scheme_factory", scheme_factories)
    test_run_problem.parametrize("problem_factory", problem_factories)
    test_run_problem.parametrize("refinement_level", params["refinements"])
    test_run_problem.parametrize("dt", params["dts"])
