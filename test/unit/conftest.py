
#from cbcflow import *
from numpy.random import random

def create_scheme_factories():

    # Make list of scheme classes that should work with default parameters:
    working_schemes = [
        IPCS,
        SegregatedIPCS,
        SegregatedIPCS_Optimized,
        IPCS_Stabilized,
        PenaltyIPCS,
        SegregatedPenaltyIPCS,
        IPCS_Stable,
        CoupledNonLinear,
        Stokes,
    #    CoupledPicard, # WIP: Needs to set pressure average
        ]

    # Print list of scheme classes that are not in the above list:
    missing_schemes = set(all_schemes) - set(scheme.__name__ for scheme in working_schemes)
    if missing_schemes:
        print
        print "Not testing schemes:"
        print '\n'.join(sorted('    '+scheme for scheme in missing_schemes))
        print

    scheme_factories = []

    # Make list of factory functions for schemes with default parameters:
    scheme_factories += [lambda: Scheme() for Scheme in working_schemes]

    # Add list of factory functions for schemes with non-default parameters:
    scheme_factories += [
        lambda: IPCS_Stable(),
        lambda: IPCS_Stabilized({'theta':0.0}),
        lambda: IPCS_Stabilized({'theta':1.0}),
        lambda: IPCS_Stabilized({'theta':0.5}),
        #lambda: IPCS_Stable({'adaptive_timestepping':True}),
        ]
    return scheme_factories

def pytest_addoption(parser):
    parser.addoption("--all", action="store_true",
        help="run all combinations")


def pytest_generate_tests(metafunc):
    if 'dim' in metafunc.fixturenames:
        metafunc.parametrize("dim", [2, 3])

    # TODO: Make options to select all or subset of schemes for this factory,
    #       copy from or look at regression conftest,
    if 'scheme_factory' in metafunc.fixturenames:
        metafunc.parametrize("scheme_factory", create_scheme_factories())
        
    if 'D' in metafunc.fixturenames:
        metafunc.parametrize("D", [2,3])
        
    
    if 'start_time' in metafunc.fixturenames:
        start_times = [0.0]
        #start_times = [0.18]
        if metafunc.config.option.all:
            start_times += list(0.8*random(3))
        metafunc.parametrize("start_time", start_times)
        
    if 'end_time' in metafunc.fixturenames:
        end_times = [2.0]
        #end_times = [1.55]
        if metafunc.config.option.all:
            end_times += list(1.2+0.8*random(3))
        metafunc.parametrize("end_time", end_times)
        
    if 'dt' in metafunc.fixturenames:
        dts = [0.1]
        #dts = [0.2]
        if metafunc.config.option.all:
            dts += [0.05+0.05*random(), 0.2+0.2*random()]
        metafunc.parametrize("dt", dts)
