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
"""A collection of postprocessing fields (Fields) to be used by a NSPostProcessor object."""

# Lists of available field names
from cbcflow.fields.basic import basic_fields
from cbcpost import Field

# Import field classes from modules with same name
for f in basic_fields:
    exec("from cbcflow.fields.basic.%s import %s" % (f, f))
#for f in meta_fields:
#    exec("from cbcflow.fields.meta.%s import %s" % (f, f))

# Make a mapping from name to type, for use in NSPostProcessor
field_classes = { f: eval(f) for f in basic_fields }
assert all(issubclass(c, Field) for c in field_classes.itervalues())

def show_fields():
    "Lists which fields are available."
    print "Postprocessing fields available by name:"
    print "\n".join("    " + f for f in basic_fields)
    print "Postprocessing fields available with parameters:"
    print "\n".join("    " + f for f in meta_fields)

