#!/usr/bin/env python
"""
Tests of the custom boundary condition functionality in cbcflow.

FIXME: This is a prime candidate for better unit tests!!!
"""

import unittest

import numpy as np
from cbcflow import make_pouseille_bcs, make_womersley_bcs
from cbcflow.dol import Function, VectorFunctionSpace, Mesh, MeshFunction, Expression, DirichletBC, assemble, ds

class TestBoundaryConditions(unittest.TestCase):
    def setUp(self):
        # Geometry
        self.mesh = Mesh("data/cylinder_4k.xml.gz")
        self.facet_domains = MeshFunction("size_t", self.mesh, 3-1, self.mesh.domains())
        self.indicator = 1

        # Temporal profile coefficients
        ts = [0.0, 0.2, 0.3, 0.4, 0.8]
        Qs = [1.0, 5.0, 3.0, 2.0, 1.0]
        self.coeffs = zip(ts,Qs)
        self.period = max(ts)

        # Function space and function
        self.V = VectorFunctionSpace(self.mesh, "CG", 1)
        self.u = Function(self.V)

    def test_womersley(self):
        nu = 4.0
        expressions = make_womersley_bcs(self.coeffs, self.mesh, self.indicator, nu, None, self.facet_domains)

        # Test that expressions are Expressions
        for bc in expressions:
            self.assertTrue(isinstance(bc, Expression))

        # Test that expressions can be applied as BCs providing nonzero change to rhs vector for some ts
        for t in np.linspace(0.0, self.period, 10):
            self.u.vector().zero()

            for bc in expressions:
                bc.set_t(t)

            for i, bc in enumerate(expressions):
                dbc = DirichletBC(self.V.sub(i), bc, self.indicator)
                dbc.apply(self.u.vector())

            self.assertGreater(self.u.vector().norm('l2'), 0.0)

        # FIXME: Now test that the values are correct!

    def test_pouseille(self):
        expressions = make_pouseille_bcs(self.coeffs, self.mesh, self.indicator, None, self.facet_domains)

        # Test that expressions are Expressions
        for bc in expressions:
            self.assertTrue(isinstance(bc, Expression))

        # Test that expressions can be applied as BCs providing nonzero change to rhs vector for some ts
        for t in np.linspace(0.0, self.period, 10):
            self.u.vector().zero()

            for bc in expressions:
                bc.set_t(t)

            for i, bc in enumerate(expressions):
                dbc = DirichletBC(self.V.sub(i), bc, self.indicator)
                dbc.apply(self.u.vector())

            self.assertGreater(self.u.vector().norm('l2'), 0.0)

        # FIXME: Now test that the values are correct!

    def test_womersley_is_pouseille_with_stationary_coefficients(self):
        ts = [0.0, 0.2, 0.3, 0.4, 0.8]
        Qs = [1.0, 1.0, 1.0, 1.0, 1.0]
        coeffs = zip(ts, Qs)
        period = max(ts)

        nu = 1.0
        wexpressions = make_womersley_bcs(coeffs, self.mesh, self.indicator, nu, None, self.facet_domains)
        pexpressions = make_pouseille_bcs(coeffs, self.mesh, self.indicator, None, self.facet_domains)

        dsi = ds[self.facet_domains](self.indicator)

        for t in np.linspace(0.0, period, 10):
            for bc in wexpressions:
                bc.set_t(t)
            for bc in pexpressions:
                bc.set_t(t)

            for wbc, pbc in zip(wexpressions, pexpressions):
                wnorm = np.sqrt(assemble(wbc**2*dsi, mesh=self.mesh))
                pnorm = np.sqrt(assemble(pbc**2*dsi, mesh=self.mesh))
                diff = np.sqrt(assemble((wbc-pbc)**2*dsi, mesh=self.mesh))
                if 0:
                    print "wnorm", wnorm
                    print "pnorm", pnorm
                    print "diff", diff
                    print "rel", diff/pnorm # |w-p|_Gamma / |p|_Gamma
                self.assertAlmostEqual(diff, 0.0)
                self.assertAlmostEqual(wnorm, pnorm)

    def test_pouseille_flow_rates_are_correct(self):
        pass # FIXME

    def test_womersley_flow_rates_are_correct(self):
        pass # FIXME
