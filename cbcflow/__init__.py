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
r"""This is the base module of cbcflow.

To use cbcflow, do::

    from cbcflow import *
"""

__version__ = "1.3.0"

# Basic utilities
from cbcflow.core.paramdict import ParamDict

from cbcflow.post import *

# Core component interfaces
from cbcflow.core.nsproblem import NSProblem
from cbcflow.core.nspostprocessor import NSPostProcessor
from cbcflow.core.nsscheme import NSScheme
from cbcflow.core.nssolver import NSSolver
from cbcflow.core.nsreplay import NSReplay

# Problem inspection utilities
from cbcflow.utils.core.show import show_problem

# Boundary condition utilities
from cbcflow.bcs import *

# Postprocessing utilities and fields
from cbcflow.fields import *

# Navier-Stokes solver schemes
from cbcflow.schemes import *

# This is too broad:
#from cbcflow.utils import *
