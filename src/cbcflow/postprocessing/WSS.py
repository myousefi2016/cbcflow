
from .PPField import PPField
from ..core.utils import cbcflow_print, cbcflow_warning
from ..core.utils import sigma

from dolfin import *

from numpy import where


def local_mesh_to_boundary_dofmap(boundary, V, Vb):
    D = boundary.topology().dim()
    mesh = V.mesh()

    V_dm = V.dofmap()
    Vb_dm = Vb.dofmap()


    mesh_tab_coords = V_dm.tabulate_all_coordinates(mesh)
    boundary_tab_coords = Vb_dm.tabulate_all_coordinates(Vb.mesh())

    dofmap_to_boundary = {}

    vertex_map = boundary.entity_map(0)
    cell_map = boundary.entity_map(D)

    for i in xrange(len(cell_map)):
        boundary_cell = Cell(boundary, i)
        mesh_facet = Facet(mesh, cell_map[i])
        mesh_cell_index = mesh_facet.entities(D+1)[0]
        mesh_cell = Cell(mesh, mesh_cell_index)

        cell_dofs = V_dm.cell_dofs(mesh_cell_index )
        boundary_dofs = Vb_dm.cell_dofs(i)

        if V_dm.num_entity_dofs(0) > 0:
            for v_idx in boundary_cell.entities(0):
                mesh_v_idx = vertex_map[int(v_idx)]

                mesh_list_idx = where(mesh_cell.entities(0) == mesh_v_idx)[0][0]
                boundary_list_idx = where(boundary_cell.entities(0) == v_idx)[0][0]

                bdofs = boundary_dofs[Vb_dm.tabulate_entity_dofs(0, boundary_list_idx)]
                cdofs = cell_dofs[V_dm.tabulate_entity_dofs(0, mesh_list_idx)]
                for bdof, cdof in zip(bdofs, cdofs):
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
        """
        FIXME: CG1 boundary gives very odd results
        degree = spaces.V.ufl_element().degree()
        if degree > 1:
            if degree > 2:
                cbcflow_warning("WSS is reduced to piecewise linears.")
            Q = VectorFunctionSpace(problem.mesh, "CG", 1)
            Q_boundary = VectorFunctionSpace(boundary, "CG", 1)
        else:
            Q = VectorFunctionSpace(problem.mesh, "DG", 0)
            Q_boundary = VectorFunctionSpace(boundary, "DG", 0)
        """

        # Create vector DG0 space on full mesh
        #Q = VectorFunctionSpace(problem.mesh, "DG", 0)
        Q = spaces.get_space(0, 1)

        # Create vector DG0 on boundary mesh
        boundary = BoundaryMesh(problem.mesh, "exterior") # TODO: Move construction to spaces?
        Q_boundary = VectorFunctionSpace(boundary, "DG", 0)

        self.v = TestFunction(Q)
        self.tau = Function(Q, name="WSS")
        self.tau_boundary = Function(Q_boundary, name="WSS")

        self.local_dofmapping = local_mesh_to_boundary_dofmap(boundary, Q, Q_boundary)
        self._keys = self.local_dofmapping.keys()
        self._values = self.local_dofmapping.values()

        self._scaling = 1 / FacetArea(problem.mesh)
        self._n = FacetNormal(problem.mesh)

    def compute(self, pp, spaces, problem):
        scaling = self._scaling
        n = self._n

        u = pp.get("Velocity")
        p = pp.get("Pressure")
        mu = Constant(problem.params.mu)

        T = -mu*dot((grad(u) + grad(u).T), n)
        Tn = dot(T, n)
        Tt = T - Tn*n

        tau_form = scaling*dot(self.v, Tt)*ds()

        assemble(tau_form, tensor=self.tau.vector())

        self.tau_boundary.vector()[self._keys] = self.tau.vector()[self._values]

        return self.tau_boundary
