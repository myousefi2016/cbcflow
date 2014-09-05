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
from cbcflow.post.spaces import SpacePool

class NSSpacePool():
    "A function space pool with custom named spaces for use with Navier-Stokes schemes."
    def __init__(self, mesh, u_degree, p_degree, u_family="auto", p_family="auto"):
        #SpacePool._init(self, mesh)
        assert isinstance(u_degree, int)
        assert isinstance(p_degree, int)
        assert isinstance(u_family, str)
        assert isinstance(p_family, str)
        
        self._spacepool = SpacePool(mesh)
        
        # FIXME: Make this explicit?
        for attr in dir(self._spacepool):
            if attr[0] != "_":
                if not hasattr(self, attr):
                    setattr(self, attr, getattr(self._spacepool, attr))
        
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

if __name__ == '__main__':
    from dolfin import *
    from cbcflow.post.spaces import get_grad_space
    
    mesh1 = UnitCubeMesh(8,8,8)
    bmesh = BoundaryMesh(mesh1, 'exterior')
    mesh2 = UnitSquareMesh(8,8)
    
    print SpacePool._existing.keys()
    
    print SpacePool(mesh1) is SpacePool(mesh1)
    print SpacePool(mesh2) is SpacePool(mesh2)
    
    pool = SpacePool(mesh1)
    print SpacePool(mesh1) is pool
    
    print SpacePool._existing.keys()
    del pool
    print SpacePool._existing.keys()
    exit()
    #spacepool = SpacePool(mesh1)
    print SpacePool._existing.keys()
    #spacepool_duplicate = SpacePool(mesh1)
    spacepool2 = SpacePool(mesh2)
    
    
    
    
    V = pool1.get_space(2,1)
    u = Function(V)
    DV = get_grad_space(u)
    
    print get_grad_space(u) is pool1.get_grad_space(V)
    exit()
    
    print SpacePool._existing.keys()
    
    nspool1 = NSSpacePool(mesh1, 2, 1)
    nspool2 = NSSpacePool(mesh1, 1, 1)
    
    print dir(nspool1)
    import ipdb; ipdb.set_trace()
    
    print SpacePool._existing.keys()
    