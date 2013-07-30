#!/usr/bin/env python
"""
Tests of the custom boundary condition functionality in headflow.

FIXME: This is a prime candidate for better unit tests!!!
"""

import unittest

from headflow import make_pouseille_bcs, make_womersley_bcs
from headflow.dol import Function, VectorFunctionSpace, Mesh, MeshFunction, Expression, DirichletBC

class TestBoundaryConditions(unittest.TestCase):
    def setUp(self):
        # Geometry
        self.mesh = Mesh("data/cylinder_4k.xml.gz")
        self.facet_domains = MeshFunction("size_t", self.mesh, 3-1, self.mesh.domains())
        self.indicator = 1

        # Temporal profile coefficients
        x = [0.0, 0.2, 0.3, 0.4, 0.8]
        y = [1.0, 5.0, 3.0, 2.0, 1.0]
        self.coeffs = zip(x,y)

        # Function space and function
        self.V = VectorFunctionSpace(self.mesh, "CG", 1)
        self.u = Function(self.V)

    def test_womersley(self):
        nu = 4.0
        expressions = make_womersley_bcs(self.coeffs, self.mesh, self.indicator, nu, None, self.facet_domains)

        # Test that expressions are Expressions and can be applied as BCs providing nonzero change to rhs vector
        self.u.vector().zero()
        for i, bc in enumerate(expressions):
            self.assertTrue(isinstance(bc, Expression))
            dbc = DirichletBC(self.V.sub(i), bc, 1)
            dbc.apply(self.u.vector())
        self.assertGreater(self.u.vector().norm('l2'), 0.0)

        # FIXME: Now test that the values are correct!

    def test_pouseille(self):
        expressions = make_pouseille_bcs(self.coeffs, self.mesh, self.indicator, None, self.facet_domains)

        # Test that expressions are Expressions and can be applied as BCs providing nonzero change to rhs vector
        self.u.vector().zero()
        for i, bc in enumerate(expressions):
            self.assertTrue(isinstance(bc, Expression))
            dbc = DirichletBC(self.V.sub(i), bc, 1)
            dbc.apply(self.u.vector())
        self.assertGreater(self.u.vector().norm('l2'), 0.0)

        # FIXME: Now test that the values are correct!
