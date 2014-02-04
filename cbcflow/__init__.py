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
__version__ = "0.1"
__date__ = "2013-11-26"
__license__  = "GNU GPL version 3 or any later version"

# Basic utilities
from cbcflow.core.paramdict import ParamDict

# Core component interfaces
from cbcflow.core.nsproblem import NSProblem
from cbcflow.core.nspostprocessor import NSPostProcessor
from cbcflow.core.nsscheme import NSScheme
from cbcflow.core.nssolver import NSSolver
from cbcflow.core.nsreplay import NSReplay

# Problem inspection utilities
from cbcflow.core.show import show_problem

# Timestepping utilities
from cbcflow.core.adaptivetimestepping import AdaptiveTimestepping
from cbcflow.core.constanttimestepping import ConstantTimestepping

# Boundary condition utilities
from cbcflow.bcs import *

# Postprocessing utilities and fields
from cbcflow.fields import *

# Navier-Stokes solver schemes
from cbcflow.schemes import *

# Import utils
#from cbcflow.utils.fenicstools import *
#from cbcflow.utils.create_slicemesh import create_slicemesh
#from cbcflow.utils.Slice import Slice
