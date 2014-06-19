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
"""These schemes are considered experimental, and one can not generally expect
correct simulations results with these.

They implement a variety of ideas for solving the Navier-Stokes equations, but
lack the testing and validation of the *official* schemes.
"""

### Experimental schemes (i.e. work in progress, not tested well, not documented)
# The simplest ipcs schemes with different optimizations
from cbcflow.schemes.experimental.ipcs_opt_seg import SegregatedIPCS_Optimized
from cbcflow.schemes.experimental.ipcs_segregated import SegregatedIPCS
from cbcflow.schemes.experimental.ipcs_stabilized import IPCS_Stabilized
from cbcflow.schemes.experimental.ipcs_abi_cn import IPCS_ABI_CN

# Schemes with support for penalty pressure BCs
from cbcflow.schemes.experimental.ipcs_penalty import PenaltyIPCS
from cbcflow.schemes.experimental.ipcs_penalty_segregated import SegregatedPenaltyIPCS

# DG schemes not working in parallell
from cbcflow.schemes.experimental.bottipietro import BottiPietro
from cbcflow.schemes.experimental.piso import PISO
from cbcflow.schemes.experimental.karper import Karper

# Coupled schemes
from cbcflow.schemes.experimental.couplednonlinear import CoupledNonLinear
from cbcflow.schemes.experimental.coupled_picard import CoupledPicard
from cbcflow.schemes.experimental.stokes import Stokes
from cbcflow.schemes.experimental.coupledpreconditioned import CoupledPreconditioned
from cbcflow.schemes.experimental.coupledpreconditioned_kam import CoupledPreconditionedKAM

from cbcflow.schemes.experimental.yosida import Yosida

# Collect all schemes in list automatically
from cbcflow.core.nsscheme import NSScheme
experimental_schemes = [k for k,v in globals().items()
                        if hasattr(v, 'mro')
                        and issubclass(v, NSScheme)
                        and v is not NSScheme]
