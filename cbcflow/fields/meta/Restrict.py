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
from cbcflow.fields.bases.MetaPPField import MetaPPField
from cbcflow.utils.common.restriction_map import restriction_map
from cbcflow.utils.common.utils import cbcflow_warning
from cbcflow.fields import PPField
from dolfin import Function, FunctionSpace, VectorFunctionSpace, TensorFunctionSpace

class Restrict(MetaPPField):
    "Restrict is used to restrict a PPField to a submesh of the mesh associated with the PPField."
    def __init__(self, field, submesh, params={}, label=None):
        PPField.__init__(self, params, label)
        
        self.submesh = submesh
        
        # Store only name, don't need the field
        if isinstance(field, PPField):
            field = field.name
        self.valuename = field
        
    @property
    def name(self):
        n = "Restriction_%s" % self.valuename
        if self.label: n += "_"+self.label
        return n
    
    def before_first_compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)
        
        if not isinstance(u, Function):
            cbcflow_warning("Do not understand how to handle datatype %s" %str(type(u)))
            return None
        
        V = u.function_space()
        element = V.ufl_element()        
        family = element.family()
        degree = element.degree()
        
        if u.rank() == 0: FS = FunctionSpace(self.submesh, family, degree)
        elif u.rank() == 1: FS = VectorFunctionSpace(self.submesh, family, degree)
        elif u.rank() == 2: FS = TensorFunctionSpace(self.submesh, family, degree)
        
        self.u = Function(FS)
            
        self.restriction_map = restriction_map(V, FS)
        
        
    
    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)
        
        if not isinstance(u, Function):
            cbcflow_warning("Do not understand how to handle datatype %s" %str(type(u)))
            return None
        
        self.u.vector()[self.restriction_map.keys()] = u.vector()[self.restriction_map.values()]
        
        return self.u
        