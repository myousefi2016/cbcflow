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

from __future__ import division

from cbcflow.dol import (FacetNormal, DirichletBC, Constant, dot, as_vector, inner,
                         MPI, mpi_comm_world, compile_extension_module, GenericMatrix,
                         GenericVector)
from numpy import array, float_
# --- Boundary condition helper functions for schemes

import numpy as np
from mpi4py import MPI as MPIpy
class DirichletBC2(DirichletBC):
    """Alternative DirichletBC class to avoid communication and recalculation of dofs"""
    def __init__(self, *args, **kwargs):
        super(DirichletBC, self).__init__(*args, **kwargs)

        V = self.function_space()
        gdim = V.mesh().geometry().dim()

        dm = V.dofmap()
        if dm.is_view():
            raise RuntimeError("Unable to handle dofmap views.")
            
        #    dm = dm.collapse(V.mesh())
        all_dofs = self.get_boundary_values().keys()

        ownership_range = dm.ownership_range()
        N = ownership_range[1]-ownership_range[0]

        dofs = [_d for _d in all_dofs if _d+ownership_range[0]<ownership_range[1]]

        # Communicate dofs to be set to owning processes
        unowned_dofs = [_d for _d in all_dofs if not _d+ownership_range[0]<ownership_range[1]]
        global_unowned_dofs = [dm.local_to_global_index(_d) for _d in unowned_dofs]

        send_to = [list()]*MPI.size(mpi_comm_world())
        for i in range(MPI.size(mpi_comm_world())):
            send_to[i] = [dm.local_to_global_index(_d) for _d in unowned_dofs if dm.off_process_owner()[_d-N] == i]

        receive = [list()]*MPI.size(mpi_comm_world())
        received = []
        for i in range(MPI.size(mpi_comm_world())):
            receive[i] = MPIpy.COMM_WORLD.scatter(send_to, root=i)
            received += receive[i]

        received_local_dofs = [_d-ownership_range[0] for _d in received]

        dofs = list(set(dofs+received_local_dofs))

        # Extract coords of dofs
        self.coords = dm.tabulate_all_coordinates(V.mesh()).reshape(N,gdim)[dofs].flatten()
        self.dofs = np.array(dofs, dtype=np.intc)
        self.global_dofs = np.array([dm.local_to_global_index(_d) for _d in self.dofs], dtype=np.intc)
        self.f = self.value()

        cpp_code = '''
        void apply_bc(const Array<int>& dofs,
                      const Array<double>& coords,
                      std::shared_ptr<const GenericFunction> f,
                      GenericVector& x)
        {
            std::size_t N = dofs.size();
            if (N==0)
                return;
            std::size_t d = f->value_size();
            std::size_t D = coords.size()/dofs.size();
            //std::cout << N << ", " << d << ", " << D << std::endl;
        
            Array<double> vals(N*d);
            Array<double> _vals(d);
            Array<double> _coords(D);
            
            //f->eval(vals, coords);
        
            for (std::size_t i=0; i<N; i++)
            {
                for (std::size_t j=0; j<D; j++)
                {
                    _coords[j] = coords[i*D+j];
                }
        
                f->eval(_vals, _coords);
                
                for (std::size_t j=0; j<d; j++)
                {
                    vals[i*d+j] = _vals[j];
                }
            }
            
            x.set_local(vals.data(), dofs.size(), dofs.data());
        
        }
        '''

        compiled_module = compile_extension_module(cpp_code)
        self._apply_bc = compiled_module.apply_bc

    def apply(self, *args):
        assert 1<=len(args)<=3

        # Call DirichletBC.apply for special cases
        if len(args)==3:
            super(DirichletBC, self).apply(*args)
        elif len(args) == 2 and isinstance(args[0], GenericVector) and isinstance(args[1], GenericVector):
            super(DirichletBC, self).apply(*args)

        if len(args) == 2:
            assert isinstance(args[0], GenericMatrix)
            assert isinstance(args[1], GenericVector)

            M = args[0]
            b = args[1]
        elif len(args) == 1:
            if isinstance(args[0], GenericMatrix):
                M = args[0]
                b = None
            elif isinstance(args[0], GenericVector):
                M = None
                b = args[0]
            else:
                raise RuntimeError()

        # Set identity on rows in matrix
        if M:            
            M.ident(self.global_dofs)

        # Get and set values in vector
        if b:
            self._apply_bc(self.dofs, self.coords, self.f, b)

def _domainargs(problem, D):
    "Helper function to pass domain args if necessary."
    if isinstance(D, int):
        return (problem.facet_domains, D)
    else:
        return (D,)

def make_velocity_bcs(problem, spaces, bcs):
    bcu_raw, bcp_raw = bcs[:2]
    bcu = [DirichletBC(spaces.V.sub(d), functions[d], *_domainargs(problem, region))
           for functions, region in bcu_raw
           for d in range(len(functions))]
    return bcu

def make_mixed_velocity_bcs(problem, spaces, bcs):
    bcu_raw, bcp_raw = bcs[:2]
    bcu = [DirichletBC(spaces.W.sub(0).sub(d), functions[d], *_domainargs(problem, region))
           for functions, region in bcu_raw
           for d in range(len(functions))]
    return bcu

def make_segregated_velocity_bcs(problem, spaces, bcs):
    bcu_raw, bcp_raw = bcs[:2]
    has_eval = True
    D = spaces.U.mesh().geometry().dim()
    x = array([0.0]*D, dtype=float_)
    v = array([0.0], dtype=float_)

    for functions, region in bcu_raw:
        for f in functions:
            try:
                f.eval(v, x)
            except:
                has_eval = False
    
    if has_eval:
        BC_class = DirichletBC2
    else:
        BC_class = DirichletBC

    bcu = [[BC_class(spaces.U, functions[d], *_domainargs(problem, region))
            for d in range(len(functions))]
           for functions, region in bcu_raw]
    return bcu

def make_pressure_bcs(problem, spaces, bcs):
    bcu_raw, bcp_raw = bcs[:2]
    has_eval = True
    D = spaces.U.mesh().geometry().dim()
    x = array([0.0]*D, dtype=float_)
    v = array([0.0], dtype=float_)

    for functions, region in bcu_raw:
        for f in functions:
            try:
                f.eval(v, x)
            except:
                has_eval = False
    
    if has_eval:
        BC_class = DirichletBC2
    else:
        BC_class = DirichletBC
    
    bcp = [BC_class(spaces.Q, function, *_domainargs(problem, region))
           for function, region in bcp_raw]
    return bcp

def make_rhs_pressure_bcs(problem, spaces, bcs, v):
    bcu_raw, bcp_raw = bcs[:2]
    ds = problem.ds
    n = FacetNormal(problem.mesh)
    Lbc = -sum(dot(function*n, v)*ds(region) for (function, region) in bcp_raw)
    return Lbc
