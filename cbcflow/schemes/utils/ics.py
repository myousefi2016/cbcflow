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

from __future__ import division

from cbcflow.dol import as_vector, project

# --- Initial condition helper functions for schemes

def assign_ics_mixed(up0, spaces, ics):
    up = as_vector(list(ics[0]) + [ics[1]])
    up0.assign(project(up, spaces.W)) #, name="up0_init_projection"))

def assign_ics_split(u0, p0, spaces, ics):
    u = as_vector(list(ics[0]))
    p = ics[1]
    u0.assign(project(u, spaces.V)) #, name="u0_init_projection"))
    p0.assign(project(p, spaces.Q)) #, name="p0_init_projection"))

def assign_ics_segregated(u0, p0, spaces, ics):
    for d in spaces.dims:
        u0[d].assign(project(ics[0][d], spaces.U)) #, name="u0_%d_init_projection"%d))
    p0.assign(project(ics[1], spaces.Q)) #, name="p0_init_projection"))
