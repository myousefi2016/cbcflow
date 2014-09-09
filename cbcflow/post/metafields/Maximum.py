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
from cbcflow.post.fieldbases.MetaField import MetaField
from dolfin import Function, MPI
import numpy

class Maximum(MetaField):
    def compute(self, get):
        u = get(self.valuename)
        
        if u == None:
            return None
        
        if isinstance(u, Function):
            return MPI.max(numpy.max(u.vector().array()))
        elif hasattr(u, "__len__"):
            return MPI.max(max(u))
        elif isinstance(u, (float,int,long)):
            return MPI.max(u)
        else:
            raise Exception("Unable to take max of %s" %str(u))
        