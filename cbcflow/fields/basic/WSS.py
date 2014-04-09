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

from cbcflow.fields.bases.PPField import PPField

from dolfin import (BoundaryMesh, VectorFunctionSpace, Cell, Facet,
                    TestFunction, Function, FacetArea, FacetNormal,
                    Constant, dot, grad, ds, assemble, inner, dx,
                    TrialFunction, Vector, solve, MPI, LinearSolver)

from numpy import where


def local_mesh_to_boundary_dofmap(boundary, V, Vb):
    "Find the mapping between dofs on boundary and boundary dofs of full mesh"
    #print "hei"
    D = boundary.topology().dim()
    mesh = V.mesh()

    V_dm = V.dofmap()
    Vb_dm = Vb.dofmap()

    dofmap_to_boundary = {}

    vertex_map = boundary.entity_map(0)
    cell_map = boundary.entity_map(D)

    for i in xrange(len(cell_map)):
        boundary_cell = Cell(boundary, i)
        mesh_facet = Facet(mesh, cell_map[i])
        mesh_cell_index = mesh_facet.entities(D+1)[0]
        mesh_cell = Cell(mesh, mesh_cell_index)

        cell_dofs = V_dm.cell_dofs(mesh_cell_index)
        boundary_dofs = Vb_dm.cell_dofs(i)

        if V_dm.num_entity_dofs(0) > 0:
            for v_idx in boundary_cell.entities(0):
                if v_idx in boundary.topology().shared_entities(0):
                    if boundary.topology().shared_entities(0)[v_idx] == MPI.process_number():
                        continue
                
                mesh_v_idx = vertex_map[int(v_idx)]

                mesh_list_idx = where(mesh_cell.entities(0) == mesh_v_idx)[0][0]
                boundary_list_idx = where(boundary_cell.entities(0) == v_idx)[0][0]

                bdofs = boundary_dofs[Vb_dm.tabulate_entity_dofs(0, boundary_list_idx)]
                cdofs = cell_dofs[V_dm.tabulate_entity_dofs(0, mesh_list_idx)]

                for bdof, cdof in zip(bdofs, cdofs):
                    if not (V_dm.ownership_range()[0] <= cdof < V_dm.ownership_range()[1]):
                        continue
                    dofmap_to_boundary[bdof] = cdof

        #
        if V_dm.num_entity_dofs(3) > 0 and V_dm.num_entity_dofs(0) == 0:
            bdofs = boundary_dofs[Vb_dm.tabulate_entity_dofs(2,0)]
            cdofs = cell_dofs[V_dm.tabulate_entity_dofs(3,0)]
            for bdof, cdof in zip(bdofs, cdofs):
                dofmap_to_boundary[bdof] = cdof

    return dofmap_to_boundary



class WSS(PPField):
    def before_first_compute(self, pp, spaces, problem):
        boundary = BoundaryMesh(problem.mesh, "exterior") # TODO: Move construction to spaces?
        degree = spaces.V.ufl_element().degree()
        if degree <= 2:
            Q = spaces.DU
        else:
            cbcflow_warning("Unable to handle higher order WSS space. Using CG1.")
            Q = spaces.get(1,1)
        Q_boundary = VectorFunctionSpace(boundary, Q.ufl_element().family(), Q.ufl_element().degree())

        self.v = TestFunction(Q)
        self.tau = Function(Q, name="WSS_full")
        self.tau_boundary = Function(Q_boundary, name="WSS")

        local_dofmapping = local_mesh_to_boundary_dofmap(boundary, Q, Q_boundary)
        self._keys = local_dofmapping.keys()
        self._values = local_dofmapping.values()
        
        Mb = assemble(inner(TestFunction(Q_boundary), TrialFunction(Q_boundary))*dx)
        self.solver = LinearSolver("gmres", "hypre_euclid")
        self.solver.set_operator(Mb)

        self.b = Function(Q_boundary).vector()
        
        self._n = FacetNormal(problem.mesh)


    def compute(self, pp, spaces, problem):
        n = self._n

        u = pp.get("Velocity")
        
        mu = Constant(problem.params.mu)

        T = -mu*dot((grad(u) + grad(u).T), n)
        Tn = dot(T, n)
        Tt = T - Tn*n
        
        tau_form = dot(self.v, Tt)*ds()
        assemble(tau_form, tensor=self.tau.vector(), reset_sparsity=False)
        
        self.b[self._keys] = self.tau.vector()[self._values]
        
        # Ensure proper scaling
        self.solver.solve(self.tau_boundary.vector(), self.b)

        return self.tau_boundary
