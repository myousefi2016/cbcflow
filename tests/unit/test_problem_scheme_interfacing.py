#!/usr/bin/env python
"""
Tests of the interaction between all the schemes
and the problem interface.

These tests should act as a documentation of
the control and data flow between the classes.
"""

import unittest

from headflow import ParamDict, NSProblem2, all_schemes
from headflow import *
from headflow.dol import *

# Make meshes only once
square = UnitSquareMesh(3,3)
square_facet_domains = FacetFunction("size_t", square)
square_facet_domains.set_all(0)

cube = UnitCubeMesh(3,3,3)
cube_facet_domains = FacetFunction("size_t", cube)
cube_facet_domains.set_all(0)

c0 = Constant(0.0)

class MockProblem2D(NSProblem2):
    def __init__(self, params=None):
        NSProblem2.__init__(self, params)
        mesh = square
        facet_domains = square_facet_domains
        self.initialize_geometry(mesh, facet_domains=facet_domains)

        self._calls = []

    @classmethod
    def default_user_params(cls):
        return ParamDict(
            T0=1.0,
            T=2.0,
            dt=0.5,

            mu=1.0,
            rho=1.0,
            )

    def observations(self, spaces, t):
        # Record this call
        self._calls.append( ("observations", type(spaces), type(t), float(t)) )

        # Return something conforming to problem specification
        return ()

    def controls(self, spaces):
        # Record this call
        self._calls.append( ("controls", type(spaces)) )

        # Return something conforming to problem specification
        return []

    def initial_conditions(self, spaces, controls):
        # Record this call
        self._calls.append( ("initial_conditions", type(spaces), len(controls)) )

        # Return something conforming to problem specification
        u0 = as_vector((c0,c0))
        p0 = c0
        ics = (u0, p0)
        return ics

    def boundary_conditions(self, spaces, u, p, t, controls):
        # Record this call
        self._calls.append( ("boundary_conditions", type(spaces), type(u), type(p), type(t), float(t), len(controls)) )

        # Return something conforming to problem specification
        d = 2
        bcu = [([c0]*d, 0)]
        bcp = [(c0, 1)]
        bcs = (bcu, bcp)
        return bcs


class TestProblemSchemeInterfacing(unittest.TestCase):
    @unittest.skip("Need a minimally computable mock problem.")
    def test_schemes_call_update_properly(self):

        schemes = [CoupledPicard]
        #schemes = all_schemes

        problems = [MockProblem2D]

        for Problem in problems:
            for Scheme in schemes:
                # Mock postprocessing update function
                update_record = []
                def update(u, p, t, timestep):
                    update_record.append((float(t), int(timestep)))

                # Run scheme with mock problem
                problem = Problem()
                scheme = Scheme()
                namespace = scheme.solve(problem, update)

                # Check that update has been called properly and that the timesteps are as they should
                self.assertEqual(len(update_record), 3)
                self.assertEqual([r[1] for r in update_record], [0,1,2])
                self.assertEqual([r[0] for r in update_record], [1.0,1.5,2.0])

                # TODO: Add checks for all problem interface components

                # TODO: Inspect problem._calls
                print [c[0] for c in problem._calls]

                # TODO: Inspect namespace contents
                self.assertIn("spaces", namespace)
                self.assertIn("observations", namespace)
