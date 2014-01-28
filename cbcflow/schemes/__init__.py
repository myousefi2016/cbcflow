"""A collection of Navier-Stokes schemes."""

# The simplest ipcs schemes with different optimizations
from .ipcs import IPCS
from .ipcs_segregated import SegregatedIPCS
from .ipcs_opt_seg import SegregatedIPCS_Optimized

# Schemes with stabilization
from .ipcs_stabilized import IPCS_Stabilized
from .ipcs_stable import IPCS_Stable

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
from .yosida import Yosida

# Collect all schemes in list automatically
from ..core.nsscheme import NSScheme
import types
all_schemes = [v for v in globals().values()
               if hasattr(v, 'mro')
               and issubclass(v, NSScheme)
               and v is not NSScheme]