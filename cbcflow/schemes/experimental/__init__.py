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

### Experimental schemes (i.e. work in progress, not tested well, not documented)
# The simplest ipcs schemes with different optimizations
from .ipcs_opt_seg import SegregatedIPCS_Optimized

# Schemes with support for penalty pressure BCs
from .ipcs_penalty import PenaltyIPCS
from .ipcs_penalty_segregated import SegregatedPenaltyIPCS

# DG schemes not working in parallell
from .bottipietro import BottiPietro
from .piso import PISO
from .karper import Karper

# Coupled schemes
from .couplednonlinear import CoupledNonLinear
from .coupled_picard import CoupledPicard
from .stokes import Stokes
from .coupledpreconditioned import CoupledPreconditoned
from .coupledpreconditioned_kam import CoupledPreconditonedKAM

# Collect all schemes in list automatically
from cbcflow.core.nsscheme import NSScheme
experimental_schemes = [k for k,v in globals().items()
                        if hasattr(v, 'mro')
                        and issubclass(v, NSScheme)
                        and v is not NSScheme]
