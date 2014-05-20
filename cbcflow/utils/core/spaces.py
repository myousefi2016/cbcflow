# Copyright (C) 2010-2014 Simula Research Laboratory
#
# This file is part of CBCFLOW.
#
# CBCFLOW is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CBCFLOW is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CBCFLOW. If not, see <http://www.gnu.org/licenses/>.


from cbcflow.dol import (FunctionSpace, VectorFunctionSpace,
                         TensorFunctionSpace, BoundaryMesh)

def galerkin_family(degree):
    return "CG" if degree > 0 else "DG"

def decide_family(family, degree):
    return galerkin_family(degree) if family == "auto" else family

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
        
        self._boundary = None

    def get_custom_space(self, family, degree, shape, boundary=False):
        if boundary:
            mesh = self.BoundaryMesh
            key = (family, degree, shape, boundary)
        else:
            mesh = self.mesh
            key = (family, degree, shape)
        space = self._spaces.get(key)
        if space is None:
            rank = len(shape)
            if rank == 0:
                space = FunctionSpace(mesh, family, degree)
            elif rank == 1:
                space = VectorFunctionSpace(mesh, family, degree, shape[0])
            else:
                space = TensorFunctionSpace(mesh, family, degree, shape)
            self._spaces[key] = space
        return space

    def get_space(self, degree, rank, family="auto", boundary=False):
        family = decide_family(family, degree)
        shape = (self.gdim,)*rank
        return self.get_custom_space(family, degree, shape, boundary)
    
    @property
    def BoundaryMesh(self):
        if self._boundary == None:
            self._boundary = BoundaryMesh(self.mesh, "exterior")
        return self._boundary
        

class NSSpacePool(SpacePool):
    "A function space pool with custom named spaces for use with Navier-Stokes schemes."
    def __init__(self, mesh, u_degree, p_degree, u_family="auto", p_family="auto"):
        SpacePool.__init__(self, mesh)
        assert isinstance(u_degree, int)
        assert isinstance(p_degree, int)
        assert isinstance(u_family, str)
        assert isinstance(p_family, str)
        self.u_degree = u_degree
        self.p_degree = p_degree
        self.u_family = u_family
        self.p_family = p_family

    @property
    def U_CG1(self): # TODO: Remove, use get_space instead
        return self.get_space(1, 0)

    @property
    def V_CG1(self): # TODO: Remove, use get_space instead
        return self.get_space(1, 1)

    @property
    def U(self):
        "Scalar valued space for velocity components."
        return self.get_space(self.u_degree, 0, family=self.u_family)

    @property
    def V(self):
        "Vector valued space for velocity vector."
        return self.get_space(self.u_degree, 1, family=self.u_family)

    @property
    def Q(self):
        "Scalar valued space for pressure."
        return self.get_space(self.p_degree, 0, family=self.p_family)

    @property
    def DU0(self):
        "Scalar valued space for gradient component of single velocity component."
        return self.get_space(self.u_degree-1, 0, family=self.u_family)

    @property
    def DU(self):
        "Vector valued space for gradients of single velocity components."
        return self.get_space(self.u_degree-1, 1, family=self.u_family)

    @property
    def DV(self):
        "Tensor valued space for gradients of velocity vector."
        return self.get_space(self.u_degree-1, 2, family=self.u_family)

    @property
    def DQ0(self):
        "Scalar valued space for pressure gradient component."
        return self.get_space(self.p_degree-1, 0, family=self.p_family)

    @property
    def DQ(self):
        "Vector valued space for pressure gradient."
        return self.get_space(self.p_degree-1, 1, family=self.p_family)

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
    @property
    def Ubc(self):
        "List of scalar valued spaces for setting velocity BCs."
        return [self.U for d in self.gdims]

    @property
    def Qbc(self):
        "Scalar valued space for setting pressure BCs."
        return self.Q
