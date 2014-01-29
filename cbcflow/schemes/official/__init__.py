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

### Schemes tested and documented well FIXME: Define criteria and make these schemes pass

# The simplest ipcs schemes with different optimizations
from .ipcs import IPCS
from .ipcs_segregated import SegregatedIPCS

# Schemes with stabilization
from .ipcs_stabilized import IPCS_Stabilized
from .ipcs_stable import IPCS_Stable

# Coupled schemes
from .yosida import Yosida

# Collect all schemes in list automatically
from ..core.nsscheme import NSScheme
import types
official_schemes = [v for v in globals().values()
                    if hasattr(v, 'mro')
                    and issubclass(v, NSScheme)
                    and v is not NSScheme]
