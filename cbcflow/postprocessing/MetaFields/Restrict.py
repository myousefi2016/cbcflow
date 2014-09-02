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
from numpy import array, uint

class Restrict(MetaPPField):
    "Restrict is used to restrict a PPField to a submesh of the mesh associated with the PPField."
    def __init__(self, field, submesh, params={}, label=None):
        MetaPPField.__init__(self, field, params, label)
        
        self.submesh = submesh
        
    @property
    def name(self):
        n = "Restrict_%s" % self.valuename
        if self.label: n += "_"+self.label
        return n
    
    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)
        
        if u == None:
            return None
        
        if not isinstance(u, Function):
            cbcflow_warning("Do not understand how to handle datatype %s" %str(type(u)))
            return None
        
        #if not hasattr(self, "restriction_map"):
        if not hasattr(self, "keys"):
            V = u.function_space()
            element = V.ufl_element()        
            family = element.family()
            degree = element.degree()
            
            if u.rank() == 0: FS = FunctionSpace(self.submesh, family, degree)
            elif u.rank() == 1: FS = VectorFunctionSpace(self.submesh, family, degree)
            elif u.rank() == 2: FS = TensorFunctionSpace(self.submesh, family, degree)
            
            self.u = Function(FS)
            
            
            #self.restriction_map = restriction_map(V, FS)
            rmap = restriction_map(V, FS)
            self.keys = array(rmap.keys(), dtype=uint)
            self.values = array(rmap.values(), dtype=uint)
            
            
        #self.u.vector()[self.restriction_map.keys()] = u.vector()[self.restriction_map.values()]
        self.u.vector()[self.keys] = u.vector()[self.values]
        return self.u
        