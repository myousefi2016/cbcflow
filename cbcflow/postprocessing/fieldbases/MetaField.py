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

from cbcflow.fields.bases.Field import Field

class MetaField(Field):
    def __init__(self, value, params=None, label=None):
        Field.__init__(self, params, label)
        self.valuename = value.name if isinstance(value, Field) else value

    @property
    def name(self):
        n = "%s_%s" % (self.__class__.__name__, self.valuename)
        if self.label: n += "_"+self.label
        return n
    
    def after_last_compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)
        if u != "N/A":
            return self.compute(pp, spaces, problem)
