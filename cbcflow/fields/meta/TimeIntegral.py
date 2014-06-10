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
from dolfin import Function
EPS = 1e-10

class TimeIntegral(MetaPPField):
    @classmethod
    def default_params(cls):
        params = MetaPPField.default_params()
        params.update(
            finalize=True,
            )
        return params
    
    def compute(self, pp, spaces, problem):
        t1 = pp.get("t")
        t0 = pp.get("t", -1)

        assert t0 <= self.params.end_time+EPS and t1 >= self.params.start_time-EPS, "Trying to compute integral outside the integration limits!"

        # Get integrand
        u1 = pp.get(self.valuename)
        u0 = pp.get(self.valuename, -1)
        
        #if u0 == "N/A" or u1 == "N/A":
        #    import ipdb; ipdb.set_trace()
        
        assert u0 != "N/A" and u1 != "N/A"
        
        # Interpolate to integration limits, if t1 or t0 is outside
        if t0 < self.params.start_time:
            if isinstance(u0, Function): start = u0.vector() + (u1.vector()-u0.vector())/(t1-t0)*(self.params.start_time-t0)
            elif hasattr(u0, "__iter__"): start = [u0[i]+(u1[i]-u0[i])/(t1-t0)*self.params.start_time-t0]
            else: start = u0 + (u1-u0)/(t1-t0)*(self.params.start_time-t0)
            t0 = self.params.start_time
        else:
            if isinstance(u0, Function): start = u0.vector()
            else: start = u0
            
        if t1 > self.params.end_time:
            if isinstance(u0, Function): end = u0.vector() + (u1.vector()-u0.vector())/(t1-t0)*(self.params.end_time-t0)
            elif hasattr(u0, "__iter__"): end = [u0[i]+(u1[i]-u0[i])/(t1-t0)*self.params.end_time-t0]
            else: end = u0 + (u1-u0)/(t1-t0)*(self.params.end_time-t0)
            t1 = self.params.end_time
        else:
            if isinstance(u1, Function): end = u1.vector()
            else: end = u1
        
        dt = t1 - t0
        if dt == 0: dt = 1e-14 # Avoid zero-division
        
        # Add to sum
        if isinstance(u0, Function):
            # Create placeholder for sum the first time
            if not hasattr(self, "_sum"):
                self._sum = Function(u0.function_space())
            # Accumulate using trapezoidal integration
            self._sum.vector().axpy(dt/2.0, start) # FIXME: Validate this, not tested!
            self._sum.vector().axpy(dt/2.0, end)
        elif hasattr(u0, "__iter__"):
            # Create placeholder for sum the first time
            if not hasattr(self, "_sum"):
                self._sum = [0.0]*len(u0)
            
            # Accumulate using trapezoidal integration
            for i in range(len(u0)):
                self._sum[i] += dt/2.0*start[i]
                self._sum[i] += dt/2.0*end[i]
        else:
            # Create placeholder for sum the first time
            if not hasattr(self, "_sum"):
                self._sum = 0.0
            # Accumulate using trapezoidal integration
            self._sum += dt/2.0*start
            self._sum += dt/2.0*end
        
        # Store bounds for sanity check
        if not hasattr(self, "T0"):
            self.T0 = t0
        self.T1 = t1
        #print "Integrated %s from %f to %f" %(self.valuename, self.T0, self.T1)
        
        if not self.params.finalize:
            return self._sum
        else:
            return None
    
    def after_last_compute(self, pp, spaces, problem):
        if not hasattr(self, "_sum"):
            return None
        
        # Integrate last timestep if integration not completed
        if self.T1 <= self.params.end_time + EPS:
            self.compute(pp, spaces, problem)
        
        #print "Integrated %s from %f to %f" %(self.valuename, self.T0, self.T1)
        #print self._sum

        return self._sum
