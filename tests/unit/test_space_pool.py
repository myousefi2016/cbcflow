
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
        if space is None
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

    @property
    def U_CG1(self): # TODO: Remove, use get_space instead
        return self.get_space(1, 0)

    @property
    def V_CG1(self): # TODO: Remove, use get_space instead
        return self.get_space(1, 1)

class NSSpacePool(SpacePool):
    "A function space pool with custom named spaces for use with Navier-Stokes schemes."
    def __init__(self, mesh, u_degree, p_degree):
        SpacePool.__init__(self, mesh)
        self.u_degree = u_degree
        self.p_degree = p_degree

    @property
    def U(self):
        "Scalar valued velocity space."
        return self.get_space(self.u_degree, 0)

    @property
    def V(self):
        "Vector valued velocity space."
        return self.get_space(self.u_degree, 1)

    @property
    def Q(self):
        "Pressure space."
        return self.get_space(self.p_degree, 0)

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

    def borrow_function(self, space, name):
        key = space
        fr = self._free[key]
        if fr:
            f = fr.pop()
            f.rename(name)
        else:
            f = Function(space, name=name)
        return f

    def return_function(self, f):
        space = f.function_space()
        f.rename("<owned by function pool>")
        key = space
        self._free[key].append(f)


import unittest

from dolfin import UnitSquareMesh
mesh = UnitSquareMesh(1,1)

class TestSpacePools(unittest.TestCase):
    def xtest_function_pool_borrow_and_return(self): # Temporarily disabled, not tested
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
            self.assertIsNot(f1, f2)
            self.assertEqual(f1.function_space().degree(), d)
            self.assertEqual(f2.function_space().degree(), d)

        # Return some of them (but hold on to the references)
        for d in (0,):
            V = sp.get_space(d, 0)
            fp.return_function(f1[d])

        # Borrow some functions matching the returned ones
        for d in (0,):
            V = sp.get_space(d, 0)
            tmp1 = fp.borrow_function(V, "xf1[%d]"%d)
            self.assertIs(f1[d], tmp)

        # Return all functions
        for d in range(n):
            V = sp.get_space(d, 0)
            fp.return_function(f1[d])
            fp.return_function(f2[d])
