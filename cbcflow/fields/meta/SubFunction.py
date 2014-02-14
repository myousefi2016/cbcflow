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
from dolfin import Function, VectorFunctionSpace, FunctionSpace, project, as_vector

def import_fenicstools():
    import cbcflow.utils.fenicstools
    return cbcflow.utils.fenicstools

class SubFunction(PPField):
    def __init__(self, field, submesh, params=None, label=None):
        PPField.__init__(self, params, label)
        
        import imp
        try:
            imp.find_module("mpi4py")
        except:
            raise ImportError("Can't find module mpi4py. This is required for SubFunction.")

        self._ft = import_fenicstools()

        self.submesh = submesh

        # Store only name, don't need the field
        if isinstance(field, PPField):
            field = field.name
        self.valuename = field

    @property
    def name(self):
        n = "SubFunction_%s" % self.valuename
        if self.label: n += "_"+self.label
        return n

    def before_first_compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)
        
        V = u.function_space()
        element = V.ufl_element()        
        family = element.family()
        degree = element.degree()
        
        if u.rank() == 1:
            FS = VectorFunctionSpace(self.submesh, family, degree)
            FS_scalar = FS.sub(0).collapse()
            self.u0 = Function(FS_scalar)
            self.u1 = Function(FS_scalar)
            self.u2 = Function(FS_scalar)

        elif u.rank() == 0:
            FS = FunctionSpace(self.submesh, family, degree)
        else:
            raise Exception("Does not support TensorFunctionSpace yet")
        
        self.u = Function(FS, name=self.name)

    def compute(self, pp, spaces, problem):       
        u = pp.get(self.valuename)

        # FIXME: This is broken for e.g. 2D, and with latest fenicstools
        
        if u.rank() == 1:
            u0, u1, u2 = u.split()
            self.u0.assign(self._ft.interpolate_nonmatching_mesh(u0, self.u0))
            self.u1.assign(self._ft.interpolate_nonmatching_mesh(u1, self.u1))
            self.u2.assign(self._ft.interpolate_nonmatching_mesh(u2, self.u2))
            self.u.assign(project(as_vector([self.u0, self.u1, self.u2])))

        elif u.rank() == 0:
            self.u.assign(self._ft.interpolate_nonmatching_mesh(u, self.u))

        return self.u
