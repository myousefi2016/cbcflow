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

from cbcflow.post.fieldbases.Field import Field

from dolfin import (TestFunction, Function,  FacetNormal,
                    Constant, dot, grad, ds, assemble, inner, dx,
                    TrialFunction, LinearSolver)

from cbcflow.utils.common.utils import cbcflow_warning
from cbcflow.utils.fields.mesh_to_boundarymesh_dofmap import mesh_to_boundarymesh_dofmap
class WSS(Field):
    def before_first_compute(self, get):
        #boundary = spaces.BoundaryMesh #BoundaryMesh(problem.mesh, "exterior") # TODO: Move construction to spaces?
        degree = spaces.V.ufl_element().degree()
        if degree <= 2:
            Q = spaces.DU
        else:
            cbcflow_warning("Unable to handle higher order WSS space. Using CG1.")
            Q = spaces.get_space(1,1)

        Q_boundary = spaces.get_space(Q.ufl_element().degree(), 1, boundary=True)

        self.v = TestFunction(Q)
        self.tau = Function(Q, name="WSS_full")
        self.tau_boundary = Function(Q_boundary, name="WSS")

        local_dofmapping = mesh_to_boundarymesh_dofmap(spaces.BoundaryMesh, Q, Q_boundary)
        self._keys = local_dofmapping.keys()
        self._values = local_dofmapping.values()
        
        Mb = assemble(inner(TestFunction(Q_boundary), TrialFunction(Q_boundary))*dx)
        #self.solver = LinearSolver("gmres", "hypre_euclid")
        self.solver = LinearSolver("gmres", "jacobi")
        self.solver.set_operator(Mb)

        self.b = Function(Q_boundary).vector()
        
        self._n = FacetNormal(problem.mesh)


    def compute(self, get):
        n = self._n

        u = get("Velocity")
        
        if isinstance(problem.params.mu, (float, int)):
            mu = Constant(problem.params.mu)
        else:
            mu = problem.params.mu

        T = -mu*dot((grad(u) + grad(u).T), n)
        Tn = dot(T, n)
        Tt = T - Tn*n
        
        tau_form = dot(self.v, Tt)*ds()
        assemble(tau_form, tensor=self.tau.vector(), reset_sparsity=False)
        
        self.b[self._keys] = self.tau.vector()[self._values]
        
        # Ensure proper scaling
        self.solver.solve(self.tau_boundary.vector(), self.b)

        return self.tau_boundary
