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

from dolfin import Function, project, Constant

class OSI(PPField):    
    def before_first_compute(self, pp, spaces, problem):
        tau = pp.get("WSS")
        
        self.osi = Function(tau.sub(0).function_space().collapse())

    def compute(self, pp, spaces, problem):
        # Requires the fields Magnitude(TimeIntegral("WSS", label="OSI")) and
        # TimeIntegral(Magnitude("WSS"), label="OSI")
        mag_ta_wss = pp.get("Magnitude_TimeIntegral_WSS_OSI")
        ta_mag_wss = pp.get("TimeIntegral_Magnitude_WSS_OSI")
        
        self.osi.assign(project(Constant(0.5)-Constant(0.5)*(mag_ta_wss/ta_mag_wss), self.osi.function_space()))
        
        return self.osi
