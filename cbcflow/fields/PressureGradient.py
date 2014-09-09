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
from dolfin import grad, Function

class PressureGradient(Field):
    def before_first_compute(self, get):
        if self.params.assemble:
            V = spaces.get_space(0, 1)
        else:
            V = spaces.DQ
        self._function = Function(V, name=self.name)

    def compute(self, get):
        p = get("Pressure")

        expr = grad(p)

        return self.expr2function(expr, self._function)