"""
Tests of the custom boundary condition functionality in headflow.
"""

import unittest

from dolfin import Function, VectorFunctionSpace, Mesh, Expression, DirichletBC
from headflow import Pouseille, Womersley

class TestBoundaryConditions(unittest.TestCase):
    def setUp(self):
        self.mesh = Mesh("cylinder_4k.xml.gz")
        x = [0.0, 0.2, 0.4, 0.6, 0.8]
        y = [1, 5, 3, 2, 1]
        self.coeffs = zip(x,y)

        self.V = VectorFunctionSpace(self.mesh, "CG", 1)
        self.u = Function(self.V)

    def _test_bcs(self, bcs):
        # TODO: Expand
        for i, bc in enumerate(bcs):
            self.assertTrue(isinstance(bc, Expression))
            dbc = DirichletBC(self.V.sub(i), bc, 1)
            dbc.apply(self.u.vector())

    def test_pouseille(self):
        bcs = Pouseille(self.coeffs, self.mesh, 1)
        self._test_bcs(bcs)

    def test_womersley(self):
        bcs = Womersley(self.coeffs, self.mesh, 1, 4.0)
        self._test_bcs(bcs)
