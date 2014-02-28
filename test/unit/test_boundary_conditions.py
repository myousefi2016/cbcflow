#!/usr/bin/env py.test
"""
Tests of the custom boundary condition functionality in cbcflow.

FIXME: This is a prime candidate for better unit tests!!!
"""

import pytest

import numpy as np

from cbcflow import make_poiseuille_bcs, make_womersley_bcs
from cbcflow.dol import Function, VectorFunctionSpace, Mesh, MeshFunction, Expression, DirichletBC, assemble, ds, SubDomain

import os
data_dir = os.path.abspath(os.path.join(os.path.split(__file__)[0], "..", "..", "cbcflow-data"))

# This used to be setUp, quick fix for refactoring to py.test:
class Data:
    def __init__(self):
        # Geometry
        self.mesh = Mesh(os.path.join(data_dir, "pipe_3k.xml.gz"))
        self.facet_domains = MeshFunction("size_t", self.mesh, 3-1, self.mesh.domains())
        self.facet_domains.set_all(0)
        self.indicator = 1
        
        class SD(SubDomain):
            def inside(self, x, on_boundary):
                return x[0] < 1e-8 and on_boundary

        SD().mark(self.facet_domains, self.indicator)
        assert max(self.facet_domains) == 1, "No domains set for facet domains"

        # Temporal profile coefficients
        ts = [0.0, 0.2, 0.3, 0.4, 0.8]
        Qs = [1.0, 5.0, 3.0, 2.0, 1.0]
        self.coeffs = zip(ts,Qs)
        self.period = max(ts)

        # Function space and function
        self.V = VectorFunctionSpace(self.mesh, "CG", 1)
        self.u = Function(self.V)

@pytest.fixture
def data():
    return Data()

def test_womersley(data):
    nu = 4.0
    expressions = make_womersley_bcs(data.coeffs, data.mesh, data.indicator, nu, None, data.facet_domains)

    # Test that expressions are Expressions
    for bc in expressions:
        assert isinstance(bc, Expression)

    # Test that expressions can be applied as BCs providing nonzero change to rhs vector for some ts
    for t in np.linspace(0.0, data.period, 10):
        data.u.vector().zero()

        for bc in expressions:
            bc.set_t(t)

        for i, bc in enumerate(expressions):
            dbc = DirichletBC(data.V.sub(i), bc, data.facet_domains, data.indicator)
            dbc.apply(data.u.vector())

        assert data.u.vector().norm('l2') > 0.0

    # FIXME: Now test that the values are correct!

def test_poiseuille(data):
    expressions = make_poiseuille_bcs(data.coeffs, data.mesh, data.indicator, None, data.facet_domains)

    # Test that expressions are Expressions
    for bc in expressions:
        assert isinstance(bc, Expression)

    # Test that expressions can be applied as BCs providing nonzero change to rhs vector for some ts
    for t in np.linspace(0.0, data.period, 10):
        data.u.vector().zero()

        for bc in expressions:
            bc.set_t(t)

        for i, bc in enumerate(expressions):
            dbc = DirichletBC(data.V.sub(i), bc, data.facet_domains, data.indicator)
            dbc.apply(data.u.vector())

        assert data.u.vector().norm('l2') > 0.0

    # FIXME: Now test that the values are correct!

def test_womersley_is_poiseuille_with_stationary_coefficients(data):
    ts = [0.0, 0.2, 0.3, 0.4, 0.8]
    Qs = [1.0, 1.0, 1.0, 1.0, 1.0]
    coeffs = zip(ts, Qs)
    period = max(ts)

    nu = 1.0
    wexpressions = make_womersley_bcs(coeffs, data.mesh, data.indicator, nu, None, data.facet_domains)
    pexpressions = make_poiseuille_bcs(coeffs, data.mesh, data.indicator, None, data.facet_domains)

    dsi = ds[data.facet_domains](data.indicator)

    for t in np.linspace(0.0, period, 10):
        for bc in wexpressions:
            bc.set_t(t)
        for bc in pexpressions:
            bc.set_t(t)

        for wbc, pbc in zip(wexpressions, pexpressions):
            wnorm = np.sqrt(assemble(wbc**2*dsi, mesh=data.mesh))
            pnorm = np.sqrt(assemble(pbc**2*dsi, mesh=data.mesh))
            diff = np.sqrt(assemble((wbc-pbc)**2*dsi, mesh=data.mesh))
            if 0:
                print "wnorm", wnorm
                print "pnorm", pnorm
                print "diff", diff
                print "rel", diff/pnorm # |w-p|_Gamma / |p|_Gamma
            assert abs(diff) < 1e-8
            assert abs(wnorm - pnorm) < 1e-8

def test_poiseuille_flow_rates_are_correct():
    pass # FIXME

def test_womersley_flow_rates_are_correct():
    pass # FIXME
