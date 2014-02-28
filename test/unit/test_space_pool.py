#!/usr/bin/env py.test

from dolfin import UnitSquareMesh
from cbcflow.utils.core import SpacePool, NSSpacePool, NSSpacePoolMixed, NSSpacePoolSplit, NSSpacePoolSegregated

def test_spacepool_base_functionality():
    mesh = UnitSquareMesh(4,4)
    p = SpacePool(mesh)
    d = 2
    spaces = []
    shapes = [(), (3,), (2,4)]
    degrees = [0,1,2]
    for shape in shapes:
        for degree in degrees:
            V = p.get_custom_space("DG", degree, shape)
            assert V.ufl_element().degree() == degree
            assert V.ufl_element().value_shape() == shape

            rank = len(shape)
            shape2 = (d,)*rank
            U = p.get_space(degree, rank)
            assert U.ufl_element().degree() == degree
            assert U.ufl_element().value_shape() == shape2

            spaces.append((V,U))

    k = 0
    for shape in shapes:
        for degree in degrees:
            V0, U0 = spaces[k]; k += 1

            V = p.get_custom_space("DG", degree, shape)
            U = p.get_space(degree, len(shape))

            assert id(V0) == id(V)
            assert id(U0) == id(U)

def test_nsspacepool_named_spaces():
    mesh = UnitSquareMesh(4,4)
    p = NSSpacePool(mesh, 2, 1)
    d = 2

    assert p.U.ufl_element().degree() == 2
    assert p.V.ufl_element().degree() == 2
    assert p.Q.ufl_element().degree() == 1
    assert p.DU.ufl_element().degree() == 1
    assert p.DV.ufl_element().degree() == 1
    assert p.DQ.ufl_element().degree() == 0
    assert p.DU0.ufl_element().degree() == 1
    assert p.DQ0.ufl_element().degree() == 0
    assert p.W.ufl_element().degree() == 2

    assert p.U.ufl_element().value_shape() == ()
    assert p.V.ufl_element().value_shape() == (d,)
    assert p.Q.ufl_element().value_shape() == ()
    assert p.DU.ufl_element().value_shape() == (d,)
    assert p.DV.ufl_element().value_shape() == (d,d)
    assert p.DQ.ufl_element().value_shape() == (d,)
    assert p.DU0.ufl_element().value_shape() == ()
    assert p.DQ0.ufl_element().value_shape() == ()
    assert p.W.ufl_element().value_shape() == (d+1,)

def test_nsspacepool_mixed_bcspaces():
    d = 2
    mesh = UnitSquareMesh(4,4)
    p = NSSpacePoolMixed(mesh, 2, 1)
    W = p.W
    Ubc = p.Ubc
    Qbc = p.Qbc
    for i in range(d):
        assert Ubc[i].component() in (i,) # subspace i of subspace 0 of mixed space
    assert Qbc.component() in (1,) # subspace 1 of mixed space

def test_nsspacepool_split_bcspaces():
    d = 2
    mesh = UnitSquareMesh(4,4)
    p = NSSpacePoolSplit(mesh, 2, 1)
    V = p.V
    Q = p.Q
    Ubc = p.Ubc
    Qbc = p.Qbc
    for i in range(d):
        assert Ubc[i].component() in (i,)
    assert id(Q) == id(Qbc)

def test_nsspacepool_segregated_bcspaces():
    d = 2
    mesh = UnitSquareMesh(4,4)
    p = NSSpacePoolSegregated(mesh, 2, 1)
    U = p.U
    Q = p.Q
    Ubc = p.Ubc
    Qbc = p.Qbc
    for i in range(d):
        assert id(Ubc[i]) == id(U)
    assert id(Q) == id(Qbc)