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
from ..bases.PPField import PPField
from dolfin import TrialFunction, TestFunction, dot, grad, DirichletBC, DomainBoundary, dx, Constant, assemble, KrylovSolver, Vector, Function, solve

class StreamFunction(PPField):
    def before_first_compute(self, pp, spaces, problem):
        u = pp.get("Velocity")
        assert len(u) == 2, "Can only compute stream function for 2D problems"
        V = spaces.U
        psi = TrialFunction(V)
        self.q = TestFunction(V)
        a = dot(grad(psi), grad(self.q))*dx()
        
        
        self.bc = DirichletBC(V, Constant(0), DomainBoundary())
        self.A = assemble(a)
        self.L = Vector()
        self.bc.apply(self.A)
        #self.solver = KrylovSolver(A, "cg")
        self.psi = Function(V)
        

    def compute(self, pp, spaces, problem):
        u = pp.get("Velocity")
        assemble(dot(u[1].dx(0)-u[0].dx(1), self.q)*dx(), tensor=self.L)
        self.bc.apply(self.L)
        
        #self.solver.solve(self.psi.vector(), self.L)
        solve(self.A, self.psi.vector(), self.L)
        
        return self.psi
        
        
        
        
        
        
     
