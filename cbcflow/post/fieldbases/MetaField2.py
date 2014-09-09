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

class MetaField2(Field):
    def __init__(self, value1, value2, params=None, name="default", label=None):
        Field.__init__(self, params, name, label)
        self.valuename1 = value1.name if isinstance(value1, Field) else value1
        self.valuename2 = value2.name if isinstance(value2, Field) else value2

    @property
    def name(self):
        n = "%s_%s_%s" % (self._name, self.valuename1, self.valuename2)
        if self.label: n += "_"+self.label
            
        return n