
# --- Function space helper functions and classes for schemes

from dolfin import FunctionSpace, VectorFunctionSpace, TensorFunctionSpace, Function

def galerkin_family(degree):
    return "CG" if degree > 0 else "DG"

class SpacePool(object):
    "A function space pool to reuse spaces across a program."
    def __init__(self, mesh):
        # Store mesh reference to create future spaces
        self.mesh = mesh

        # Get dimensions for convenience
        cell = mesh.ufl_cell()
        self.gdim = cell.geometric_dimension()
        self.tdim = cell.topological_dimension()
        self.gdims = range(self.gdim)
        self.tdims = range(self.tdim)

        # For compatibility, remove when code has been converted
        self.d = self.gdim
        self.dims = self.gdims

        # Start with empty cache
        self._spaces = {}

    def get_custom_space(self, family, degree, shape):
        key = (family, degree, shape)
        space = self._spaces.get(key)
        if space is None:
            rank = len(shape)
            if rank == 0:
                space = FunctionSpace(self.mesh, family, degree)
            elif rank == 1:
                space = VectorFunctionSpace(self.mesh, family, degree, shape[0])
            else:
                space = TensorFunctionSpace(self.mesh, family, degree, shape)
            self._spaces[key] = space
        return space

    def get_space(self, degree, rank):
        family = galerkin_family(degree)
        shape = (self.gdim,)*rank
        return self.get_custom_space(family, degree, shape)

class NSSpacePool(SpacePool):
    "A function space pool with custom named spaces for use with Navier-Stokes schemes."
    def __init__(self, mesh, u_degree, p_degree):
        SpacePool.__init__(self, mesh)
        self.u_degree = u_degree
        self.p_degree = p_degree

    @property
    def U_CG1(self): # TODO: Remove, use get_space instead
        return self.get_space(1, 0)

    @property
    def V_CG1(self): # TODO: Remove, use get_space instead
        return self.get_space(1, 1)

    @property
    def U(self):
        "Scalar valued space for velocity components."
        return self.get_space(self.u_degree, 0)

    @property
    def V(self):
        "Vector valued space for velocity vector."
        return self.get_space(self.u_degree, 1)

    @property
    def Q(self):
        "Scalar valued space for pressure."
        return self.get_space(self.p_degree, 0)

    @property
    def DU0(self):
        "Scalar valued space for gradient component of single velocity component."
        return self.get_space(self.u_degree-1, 0)

    @property
    def DU(self):
        "Vector valued space for gradients of single velocity components."
        return self.get_space(self.u_degree-1, 1)

    @property
    def DV(self):
        "Tensor valued space for gradients of velocity vector."
        return self.get_space(self.u_degree-1, 2)

    @property
    def DQ0(self):
        "Scalar valued space for pressure gradient component."
        return self.get_space(self.p_degree-1, 0)

    @property
    def DQ(self):
        "Vector valued space for pressure gradient."
        return self.get_space(self.p_degree-1, 1)

    @property
    def W(self):
        "Mixed velocity-pressure space."
        space = self._spaces.get("W")
        if space is None:
            space = self.V*self.Q
            self._spaces["W"] = space
        return space

class NSSpacePoolMixed(NSSpacePool):
    "A function space pool with custom named spaces for use with mixed Navier-Stokes schemes."
    def __init__(self, mesh, u_degree, p_degree):
        NSSpacePool.__init__(self, mesh, u_degree, p_degree)

    @property
    def Ubc(self):
        "List of scalar valued spaces for setting velocity BCs."
        return [self.W.sub(0).sub(d) for d in self.gdims]

    @property
    def Qbc(self):
        "Scalar valued space for setting pressure BCs."
        return self.W.sub(1)

class NSSpacePoolSplit(NSSpacePool):
    "A function space pool with custom named spaces for use with split Navier-Stokes schemes."
    def __init__(self, mesh, u_degree, p_degree):
        NSSpacePool.__init__(self, mesh, u_degree, p_degree)

    @property
    def Ubc(self):
        "List of scalar valued spaces for setting velocity BCs."
        return [self.V.sub(d) for d in self.gdims]

    @property
    def Qbc(self):
        "Scalar valued space for setting pressure BCs."
        return self.Q

class NSSpacePoolSegregated(NSSpacePool):
    "A function space pool with custom named spaces for use with segregated Navier-Stokes schemes."
    def __init__(self, mesh, u_degree, p_degree):
        NSSpacePool.__init__(self, mesh, u_degree, p_degree)

    @property
    def Ubc(self):
        "List of scalar valued spaces for setting velocity BCs."
        return [self.U for d in self.gdims]

    @property
    def Qbc(self):
        "Scalar valued space for setting pressure BCs."
        return self.Q

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
