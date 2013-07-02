#!/usr/bin/env python
"""
Tests of the interaction between all the schemes
and the problem interface.

These tests should act as a documentation of
the control and data flow between the classes.
"""

import unittest

from headflow import ParamDict, NSProblem, all_schemes
from headflow import *
from headflow.dol import *

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] <= DOLFIN_EPS

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] >= 1.0-DOLFIN_EPS

# Make meshes only once
square = UnitSquareMesh(8,8)
square_facet_domains = FacetFunction("size_t", square)
square_facet_domains.set_all(0)
Left().mark(square_facet_domains, 1)
Right().mark(square_facet_domains, 2)

cube = UnitCubeMesh(5,5,5)
cube_facet_domains = FacetFunction("size_t", cube)
cube_facet_domains.set_all(0)
Left().mark(cube_facet_domains, 1)
Right().mark(cube_facet_domains, 2)

c0 = Constant(0.0)
c1 = Constant(0.0)

class MockProblem(NSProblem):
    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Select mesh and domains based on dimension parameter
        d = self.params.d
        if d == 2:
            mesh = square
            facet_domains = square_facet_domains
        elif d == 3:
            mesh = cube
            facet_domains = cube_facet_domains

        # Required initialization
        self.initialize_geometry(mesh, facet_domains=facet_domains)

        # These will be returned from problem interface functions
        self._observations = []
        self._controls = []
        self._ics = (as_vector([c0]*d),
                     c0)
        self._bcs = ([([c0]*d, 0)],
                     [(c0, 1), (c0, 2)])

        # List of recorded calls through the problem interface
        self._calls = []

    @classmethod
    def default_user_params(cls):
        return ParamDict(
            T0=1.0,
            T=2.0,
            dt=0.5,

            mu=1.0,
            rho=1.0,

            d=2,
            )

    def observations(self, spaces, t):
        # Record this call
        self._calls.append( ("observations", float(t)) )

        # Check that input data conforms to expected interface
        assert hasattr(spaces, 'U')
        assert isinstance(t, Constant)

        # Return something conforming to problem specification
        return self._observations

    def controls(self, spaces):
        # Record this call
        self._calls.append( ("controls",) )

        # Check that input data conforms to expected interface
        assert hasattr(spaces, 'V')

        # Return something conforming to problem specification
        return self._controls

    def initial_conditions(self, spaces, controls):
        # Record this call
        self._calls.append( ("initial_conditions",) )

        # Check that input data conforms to expected interface
        d = self.params.d
        assert controls is self._controls
        assert spaces.d == d

        # Return something conforming to problem specification
        return self._ics

    def boundary_conditions(self, spaces, u, p, t, controls):
        # Record this call
        self._calls.append( ("boundary_conditions", type(spaces), type(u), type(p), type(t), float(t), len(controls)) )

        # Check that input data conforms to expected interface
        assert controls is self._controls
        assert hasattr(spaces, 'Q')

        # TODO: Checks ...

        # Return something conforming to problem specification
        return self._bcs


class TestProblemSchemeInterfacing(unittest.TestCase):
    #@unittest.skip("Need a minimally computable mock problem.")
    def test_schemes_call_update_properly(self):

        # TODO: Goal:
        #schemes = all_schemes

        
        #schemes = [CoupledPicard]

        # Currently passing:
        schemes = []
        schemes.append(IPCS)
        schemes.append(SegregatedIPCS)
        schemes.append(SegregatedIPCS_Optimized)
        schemes.append(IPCS_Stabilized)
        schemes.append(PenaltyIPCS)
        schemes.append(SegregatedPenaltyIPCS)

        #schemes.append(IPCS_Stable) # Not passing, timestepping thingy
        #schemes.append(CoupledNonLinear)
        #schemes.append(CoupledPicard) # WIP: Needs to set pressure average
        #schemes.append(Stokes)

        missing_schemes = set(all_schemes) - set(schemes)
        if missing_schemes:
            print
            print "Not testing schemes:"
            print '\n'.join(sorted('    '+scheme.__name__ for scheme in missing_schemes))
            print

        problems = [MockProblem]

        for d in [2,3]:
            for Problem in problems:
                for Scheme in schemes:

                    # Mock postprocessing update function
                    update_record = []
                    def update(u, p, t, timestep, spaces):
                        update_record.append((float(t), int(timestep)))

                    # Run scheme with mock problem
                    problem = Problem({'d':d})
                    scheme = Scheme()
                    namespace = scheme.solve(problem, update)

                    # Check that update has been called properly and that the timesteps are as they should
                    self.assertEqual(len(update_record), 3)
                    self.assertEqual([r[1] for r in update_record], [0,1,2])
                    self.assertEqual([r[0] for r in update_record], [1.0,1.5,2.0])

                    # TODO: Add checks for all problem interface components

                    # Check that all problem interface functions were called
                    callnames = [c[0] for c in problem._calls]
                    self.assertIn("observations", callnames)
                    self.assertIn("controls", callnames)
                    self.assertIn("initial_conditions", callnames)
                    self.assertIn("boundary_conditions", callnames)

                    # TODO: Inspect problem._calls data

                    # Check that the returned namespace contains all expected values
                    self.assertIn("spaces", namespace)
                    self.assertIn("observations", namespace)
                    self.assertIn("controls", namespace)
                    self.assertIn("states", namespace)
                    self.assertIn("t", namespace)
                    self.assertIn("timesteps", namespace)

                    # TODO: Inspect namespace contents

                    # Check that the spaces object has the right function space properties
                    spaces = namespace["spaces"]
                    self.assertEqual(spaces.V.ufl_element().degree(), scheme.params.u_degree)
                    self.assertEqual(spaces.Q.ufl_element().degree(), scheme.params.p_degree)
                    self.assertEqual(spaces.V.ufl_element().value_shape(), (d,))
                    self.assertEqual(spaces.Q.ufl_element().value_shape(), ())
