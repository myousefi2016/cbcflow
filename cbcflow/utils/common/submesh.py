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
from cbcflow.utils.common.mpi_utils import (broadcast, distribute_meshdata,
                                            distribution, gather)
from dolfin import MPI, Mesh, MeshEditor
import numpy as np  

def create_submesh(mesh, markers, marker):
    "This function allows for a SubMesh-equivalent to be created in parallel"

    base_cell_indices = np.where(markers.array() == marker)[0]
    base_cells = mesh.cells()[base_cell_indices]
    base_vertex_indices = np.unique(base_cells.flatten())
    
    base_global_vertex_indices = sorted([mesh.topology().global_indices(0)[vi] for vi in base_vertex_indices])
    
    gi = mesh.topology().global_indices(0)
    shared_global_indices = set(base_vertex_indices).intersection(set(mesh.topology().shared_entities(0).keys()))
    shared_global_indices = [gi[vi] for vi in shared_global_indices]
    
    unshared_global_indices = list(set(base_global_vertex_indices)-set(shared_global_indices))
    unshared_vertices_dist = distribution(len(unshared_global_indices))
    
    # Number unshared vertices on separate process
    idx = sum(unshared_vertices_dist[:MPI.process_number()])
    base_to_sub_global_indices = {}
    for gi in unshared_global_indices:
        base_to_sub_global_indices[gi] = idx
        idx += 1
    
    # Gather all shared process on process 0 and assign global index
    all_shared_global_indices = gather(shared_global_indices, on_process=0, flatten=True)
    all_shared_global_indices = np.unique(all_shared_global_indices)

    shared_base_to_sub_global_indices = {}
    idx = int(MPI.max(float(max(base_to_sub_global_indices.values()+[-1e16])))+1)
    if MPI.process_number() == 0:
        for gi in all_shared_global_indices:
            shared_base_to_sub_global_indices[int(gi)] = idx
            idx += 1
    
    # Broadcast global numbering of all shared vertices    
    shared_base_to_sub_global_indices = dict(zip(broadcast(shared_base_to_sub_global_indices.keys(), 0),
                                                broadcast(shared_base_to_sub_global_indices.values(), 0)))
    
    # Join shared and unshared numbering in one dict
    base_to_sub_global_indices = dict(base_to_sub_global_indices.items()+
                                     shared_base_to_sub_global_indices.items())
    
    # Create mapping of local indices
    base_to_sub_local_indices = dict(zip(base_vertex_indices, range(len(base_vertex_indices))))

    # Define sub-cells
    sub_cells = [None]*len(base_cells)
    for i, c in enumerate(base_cells):
        sub_cells[i] = [base_to_sub_local_indices[j] for j in c]

    # Store vertices as sub_vertices[local_index] = (global_index, coordinates)
    sub_vertices = {}
    for base_local, sub_local in base_to_sub_local_indices.items():
        sub_vertices[sub_local] = (base_to_sub_global_indices[mesh.topology().global_indices(0)[base_local]],
                               mesh.coordinates()[base_local])
    
    ## Done with base mesh

    # Distribute meshdata on (if any) empty processes
    sub_cells, sub_vertices = distribute_meshdata(sub_cells, sub_vertices)
    global_cell_distribution = distribution(len(sub_cells))
    global_vertex_distribution = distribution(len(sub_vertices))

    global_num_cells = MPI.sum(len(sub_cells))
    global_num_vertices = int(MPI.max(max(base_to_sub_global_indices.values())))+1

    # Build mesh
    submesh = Mesh()
    mesh_editor = MeshEditor()
    mesh_editor.open(submesh,
                     mesh.ufl_cell().cellname(),
                     mesh.ufl_cell().topological_dimension(),
                     mesh.ufl_cell().geometric_dimension())
    
    mesh_editor.init_vertices(len(sub_vertices))
    mesh_editor.init_cells(len(sub_cells))
    
    for index, cell in enumerate(sub_cells):
        mesh_editor.add_cell(index, *cell)

    for local_index, (global_index, coordinates) in sub_vertices.items():
        #print coordinates
        mesh_editor.add_vertex_global(int(local_index), int(global_index), coordinates)

    mesh_editor.close()

    submesh.topology().init_global(0, global_num_vertices)
    submesh.topology().init_global(mesh.ufl_cell().topological_dimension(), global_num_cells)

    return submesh
        
        
if __name__ == '__main__':
    from dolfin import (UnitCubeMesh, BoundaryMesh, MeshFunction, FunctionSpace, SubMesh,
                        Expression, project, File, SubDomain, dx, assemble, Constant)
    #mesh = UnitCubeMesh(3,1,1)
    N = 16
    mesh = UnitCubeMesh(N,N,N)
    #print mesh.num_cells()
    #exit()
    mesh = BoundaryMesh(mesh, "exterior")
    #mesh = UnitSquareMesh(4,4)
    mf = MeshFunction("size_t", mesh, mesh.ufl_cell().topological_dimension())
    mf.set_all(0)
    
    class Left(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] < 0.4
        
    Left().mark(mf, 1)
    #plot(mf)
    #interactive()
    #tic()
    
    
    #tic()
    if MPI.num_processes() == 1:
        submesh = SubMesh(mesh, mf, 1)
    else:
        submesh = create_submesh(mesh, mf, 1)
    #print "toc: ", toc()
    #exit()
    #print submesh
    #plot(submesh)
    #interactive()
    
    #print "Num vertices: ", submesh.size_global(0)
    #print "Num cells: ", submesh.size_global(3)
    
    #print submesh.topology().shared_entities(0)
    
    #File("submesh.xdmf") << submesh
    
    
    V = FunctionSpace(submesh, "CG", 1)
    expr = Expression("x[0]*x[1]*x[1]+4*x[2]")
    u = project(expr, V)
    
    MPI.barrier()
    
    
    #print "hei: ", assemble(u*dx), "(0.685185185185)"
    #print "Num vertices: ", submesh.size_global(0)
    #print "Num cells: ", submesh.size_global(3)
    
    s0 = submesh.size_global(0)
    s3 = submesh.size_global(submesh.ufl_cell().topological_dimension())
    a = assemble(u*dx)
    v = assemble(Constant(1)*dx, mesh=submesh)
    if MPI.process_number() == 0:
        print "Num vertices: ", s0
        print "Num cells: ", s3
        print "assemble(u*dx): ", a
        print "Volume: ", v
    #u = Function(V)
    File("u.pvd") << u
    
    
    