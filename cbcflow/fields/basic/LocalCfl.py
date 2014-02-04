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
from dolfin import assemble, dx, sqrt, TestFunction, Function, Constant

class LocalCfl(PPField):
    def before_first_compute(self, pp, spaces, problem):
        DG0 = spaces.get_space(0, 0)
        self._v = TestFunction(DG0)
        self._cfl = Function(DG0)

    def compute(self, pp, spaces, problem):
        t1 = pp.get("t")
        t0 = pp.get("t", -1)
        dt = Constant(t1 - t0)
        u = pp.get("Velocity")

        cell = problem.mesh.ufl_cell()
        hF = cell.circumradius
        hK = cell.volume
        scaling = 1.0 / hK

        assemble((dt * sqrt(u**2) / hF)*self._v*scaling*dx(),
                 tensor=self._cfl.vector())

        return self._cfl
