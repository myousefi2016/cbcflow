
from cbcflow.core.spaces import SpacePool, NSSpacePool, NSSpacePoolMixed, NSSpacePoolSplit, NSSpacePoolSegregated
from dolfin import FunctionSpace, VectorFunctionSpace, TensorFunctionSpace, Function

from collections import defaultdict

class FunctionPool:
    def __init__(self):
        self._free = defaultdict(list)

    def _key(self, space):
        return (space.ufl_element(), space.mesh().id())

    def borrow_function(self, space, name):
        key = self._key(space)
        fr = self._free[key]
        if fr:
            f = fr.pop()
            f.rename(name, "<function borrowed from function pool>")
        else:
            f = Function(space, name=name)
        return f

    def return_function(self, f):
        space = f.function_space()
        f.rename("unused", "<function owned by function pool>")
        key = self._key(space)
        self._free[key].append(f)

import unittest

from dolfin import UnitSquareMesh
mesh = UnitSquareMesh(1,1)

class TestSpacePools(unittest.TestCase):
    def test_spacepool_base_functionality(self):
        p = SpacePool(mesh)
        d = 2
        spaces = []
        shapes = [(), (3,), (2,4)]
        degrees = [0,1,2]
        for shape in shapes:
            for degree in degrees:
                V = p.get_custom_space("DG", degree, shape)
                self.assertEqual(V.ufl_element().degree(), degree)
                self.assertEqual(V.ufl_element().value_shape(), shape)

                rank = len(shape)
                shape2 = (d,)*rank
                U = p.get_space(degree, rank)
                self.assertEqual(U.ufl_element().degree(), degree)
                self.assertEqual(U.ufl_element().value_shape(), shape2)

                spaces.append((V,U))

        k = 0
        for shape in shapes:
            for degree in degrees:
                V0, U0 = spaces[k]; k += 1

                V = p.get_custom_space("DG", degree, shape)
                U = p.get_space(degree, len(shape))

                self.assertEqual(id(V0), id(V))
                self.assertEqual(id(U0), id(U))

    def test_nsspacepool_named_spaces(self):
        p = NSSpacePool(mesh, 2, 1)
        d = 2

        self.assertEqual(p.U.ufl_element().degree(), 2)
        self.assertEqual(p.V.ufl_element().degree(), 2)
        self.assertEqual(p.Q.ufl_element().degree(), 1)
        self.assertEqual(p.DU.ufl_element().degree(), 1)
        self.assertEqual(p.DV.ufl_element().degree(), 1)
        self.assertEqual(p.DQ.ufl_element().degree(), 0)
        self.assertEqual(p.DU0.ufl_element().degree(), 1)
        self.assertEqual(p.DQ0.ufl_element().degree(), 0)
        self.assertEqual(p.W.ufl_element().degree(), 2)

        self.assertEqual(p.U.ufl_element().value_shape(), ())
        self.assertEqual(p.V.ufl_element().value_shape(), (d,))
        self.assertEqual(p.Q.ufl_element().value_shape(), ())
        self.assertEqual(p.DU.ufl_element().value_shape(), (d,))
        self.assertEqual(p.DV.ufl_element().value_shape(), (d,d))
        self.assertEqual(p.DQ.ufl_element().value_shape(), (d,))
        self.assertEqual(p.DU0.ufl_element().value_shape(), ())
        self.assertEqual(p.DQ0.ufl_element().value_shape(), ())
        self.assertEqual(p.W.ufl_element().value_shape(), (d+1,))

    def test_nsspacepool_mixed_bcspaces(self):
        d = 2
        p = NSSpacePoolMixed(mesh, 2, 1)
        W = p.W
        Ubc = p.Ubc
        Qbc = p.Qbc
        for i in range(d):
            self.assertEqual(Ubc[i].component(), (i,)) # subspace i of subspace 0 of mixed space
        self.assertEqual(Qbc.component(), (1,)) # subspace 1 of mixed space

    def test_nsspacepool_split_bcspaces(self):
        d = 2
        p = NSSpacePoolSplit(mesh, 2, 1)
        V = p.V
        Q = p.Q
        Ubc = p.Ubc
        Qbc = p.Qbc
        for i in range(d):
            self.assertEqual(Ubc[i].component(), (i,))
        self.assertEqual(id(Q), id(Qbc))

    def test_nsspacepool_segregated_bcspaces(self):
        d = 2
        p = NSSpacePoolSegregated(mesh, 2, 1)
        U = p.U
        Q = p.Q
        Ubc = p.Ubc
        Qbc = p.Qbc
        for i in range(d):
            self.assertEqual(id(Ubc[i]), id(U))
        self.assertEqual(id(Q), id(Qbc))

    def test_function_pool_borrow_and_return(self):
        # Setup pools to test
        sp = SpacePool(mesh)
        fp = FunctionPool()

        # Some test variables
        f1 = {}
        f2 = {}
        n = 3

        # Borrow some functions
        for d in range(n):
            V = sp.get_space(d, 0)
            f1[d] = fp.borrow_function(V, "f1[%d]"%d)
            f2[d] = fp.borrow_function(V, "f2[%d]"%d)
            self.assertIsNot(f1[d], f2[d])
            self.assertEqual(f1[d].element().degree(), d)
            self.assertEqual(f2[d].element().degree(), d)

        # Return some of them (but hold on to the references)
        for d in (0,):
            V = sp.get_space(d, 0)
            fp.return_function(f1[d])

        # Borrow some functions matching the returned ones
        for d in (0,):
            V = sp.get_space(d, 0)
            tmp = fp.borrow_function(V, "xf1[%d]"%d)
            self.assertIs(f1[d], tmp)

        # Return all functions
        for d in range(n):
            V = sp.get_space(d, 0)
            fp.return_function(f1[d])
            fp.return_function(f2[d])
