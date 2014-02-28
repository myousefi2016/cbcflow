
from cbcflow import *

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


def pytest_generate_tests(metafunc):
    if 'dim' in metafunc.fixturenames:
        metafunc.parametrize("dim", [2, 3])

    # TODO: Make options to select all or subset of schemes for this factory,
    #       copy from or look at regression conftest,
    if 'scheme_factory' in metafunc.fixturenames:
        metafunc.parametrize("scheme_factory", create_scheme_factories())
