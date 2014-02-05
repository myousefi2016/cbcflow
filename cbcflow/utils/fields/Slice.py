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
from dolfin import *
import numpy as np

cpp_code = '''
namespace dolfin {
    std::vector<double> distribute_vertices(int from_process, double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2)
    {
        int this_process = dolfin::MPI::process_number();
        
        static std::vector<double> points(9);
        
        if(this_process == from_process) {
            points[0] = x0;
            points[1] = y0;
            points[2] = z0;
            points[3] = x1;
            points[4] = y1;
            points[5] = z1;
            points[6] = x2;
            points[7] = y2;
            points[8] = z2;
        }    
        
        dolfin::MPI::broadcast(points, from_process);
        return points;
    }
    
    std::vector<unsigned int> distribution(int number)
    {
        // Variables to help in synchronization
        int num_processes = dolfin::MPI::num_processes();
        int this_process = dolfin::MPI::process_number();
        
        static std::vector<uint> distribution(num_processes);
    
        for(uint i=0; i<num_processes; i++) {
            if(i==this_process) {
                distribution[i] = number;
            }
            dolfin::MPI::barrier();
            dolfin::MPI::broadcast(distribution, i);    
        }
        return distribution;
  }
}
'''

cpp_module = compile_extension_module(cpp_code, additional_system_headers=["dolfin/common/MPI.h"])

class Slice(Mesh):
    def __init__(self, basemesh, point, normal):
        Mesh.__init__(self)
        
        P = np.array([point[0], point[1], point[2]])
        self.P = Constant((P[0], P[1],P[2]))

        # Create unit normal
        n = np.array([normal[0],normal[1], normal[2]])
        n = n/np.linalg.norm(n)
        self.n = Constant((n[0], n[1], n[2]))

        # Calculate the distribution of vertices around the plane
        # (sign of np.dot(p-P, n) determines which side of the plane p is on)
        vsplit = np.dot(basemesh.coordinates()-P, n)

        # Count each cells number of vertices on the "positive" side of the plane
        # Only cells with vertices on both sides of the plane intersect the plane
        operator = np.less
        npos = np.sum(vsplit[basemesh.cells()] < 0, 1)
        intersection_cells = basemesh.cells()[(npos > 0) & (npos < 4)]
        
        if len(intersection_cells) == 0:
            # Try to put "zeros" on other side of plane
            # FIXME: handle cells with vertices exactly intersecting the plane in a more robust manner.
            operator = np.greater
            npos = np.sum(vsplit[basemesh.cells()] > 0, 1)
            intersection_cells = basemesh.cells()[(npos > 0) & (npos < 4)]

        def add_vertex(vertices, p):
            p = np.atleast_2d(p)
            vertices = np.append(vertices, p, 0)
            return vertices, len(vertices)-1
            
        def add_cell(cells, cell):
            # Split cell into triangles
            for i in xrange(len(cell)-2):
                cells.append(cell[i:i+3])
        
        vertices = np.zeros((0, 3))
        cells = []
        indexes = {}
        for c in intersection_cells:
            a = operator(vsplit[c], 0)
            positives = c[np.where(a==True)[0]]
            negatives = c[np.where(a==False)[0]]
            
            cell = []
            for pp_ind in positives:
                pp = basemesh.coordinates()[pp_ind]

                for pn_ind in negatives:
                    pn = basemesh.coordinates()[pn_ind]
                    if (pp_ind, pn_ind) in indexes:
                        index = indexes[(pp_ind, pn_ind)]
                    else:
                        # Calculate intersection point with the plane
                        d = np.dot(P-pp, n)/np.dot(pp-pn, n)
                        ip = pp+(pp-pn)*d
                            
                        vertices, index = add_vertex(vertices, ip)
                        indexes[(pp_ind, pn_ind)] = index
                    
                    cell.append(index)
                    
            add_cell(cells, cell)
        MPI.barrier()

        global_num_cells = MPI.sum(len(cells))
        global_num_vertices = MPI.sum(len(vertices))
        
        global_cell_distribution = cpp_module.distribution(len(cells))
        global_vertex_distribution = cpp_module.distribution(len(vertices))
        
        # Return empty mesh if no intersections were found
        if global_num_cells == 0:
            mesh_editor = MeshEditor()
            mesh_editor.open(self, "triangle", 2, 3)
            
            mesh_editor.init_vertices(0)
            mesh_editor.init_cells(0)
            
            mesh_editor.close()
            return

        # Because dolfin does not support a distributed mesh that is empty on some processes,
        # we move a single cell from the process with the largest mesh to all processes with
        # empty meshes.
        while 0 in global_cell_distribution:
            to_process = list(global_cell_distribution).index(0)
            from_process = list(global_cell_distribution).index(max(global_cell_distribution))
            v_out = (0,)*9       
            if MPI.process_number() == from_process:
                v_out = vertices[cells[0],0:3].flatten()
            
            v_in = cpp_module.distribute_vertices(from_process, *v_out)
            
            if MPI.process_number() == to_process:
                cell = []
                for p in v_in.reshape(3,3):
                    vertices, index = add_vertex(vertices, p)
                    cell.append(index)
                cells.append(cell)
            if MPI.process_number() == from_process:
                # Remove vertices no longer used in remaining cells.
                for i in xrange(3):
                    v = cells[0][i]
                    if not any([v in c for c in cells[1:]]):
                        vertices = np.delete(vertices, v, axis=0)                        
                        for i in xrange(len(cells)):
                            cells[i] = [vi-1 if vi > v else vi for vi in cells[i]]
                cells.pop(0)
            MPI.barrier()
            # Update
            global_cell_distribution = cpp_module.distribution(len(cells))
            global_vertex_distribution = cpp_module.distribution(len(vertices))
        
        # Build mesh
        mesh_editor = MeshEditor()
        mesh_editor.open(self, "triangle", 2, 3)
        
        mesh_editor.init_vertices(len(vertices))
        mesh_editor.init_cells(len(cells))
        
        i = MPI.process_number()
        global_vertex_indices = range(sum(global_vertex_distribution[:i]), sum(global_vertex_distribution[:i+1]))
        
        for index, cell in enumerate(cells):
            mesh_editor.add_cell(index, cell[0], cell[1], cell[2])
        
        for index, vertex in enumerate(vertices):
            mesh_editor.add_vertex_global(index, global_vertex_indices[index], vertex)
        
        mesh_editor.close()
        
        self.topology().init_global(0, global_num_vertices)
        self.topology().init_global(2, global_num_cells)
        
        
if __name__ == '__main__':
    mesh = UnitCubeMesh(8,8,8)
    p = np.array([0.5, 0.3, 0.7])
    n = np.array([3,2,1])
    
    slicemesh = Slice(mesh, p, n)
    
    print slicemesh
    plot(slicemesh, interactive=True)
    