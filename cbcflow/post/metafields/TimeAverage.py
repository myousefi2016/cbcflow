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
#from cbcflow.fields.meta.TimeIntegral import TimeIntegral
from cbcflow.post.metafields.TimeIntegral import TimeIntegral
from dolfin import Function

class TimeAverage(TimeIntegral):
    def compute(self, get):
        
        ti = super(TimeAverage, self).compute(get)
        
        if self.params.finalize:
            return None
        else:
            return self.scale(ti)
        
        # Make sure dependencies are read by postprocessor
        # (This code is never reached, just inspected)
        get("t")
        get("t", -1)
        get(self.valuename)
        get(self.valuename, -1)
        
    def after_last_compute(self, get):
        ti = super(TimeAverage, self).after_last_compute(get)
        
        ta = self.scale(ti)
        #print "Averaged %s from %f (start_time=%f) to %f (end_time=%f) (result=%s)" %(self.valuename, self.T0, self.params.start_time, self.T1, self.params.end_time,  str(ta))
        
        return ta


    def scale(self, ti):
        if ti == None:
            return None
        
        scale_factor = 1.0/(self.T1-self.T0)
        
        if isinstance(ti, Function):
            ta = Function(ti)
            ta.vector()[:] *= scale_factor
        elif hasattr(ti, "__len__"):
            ta = [scale_factor*_ti for _ti in ti]
        else:
            ta = scale_factor*ti

        return ta