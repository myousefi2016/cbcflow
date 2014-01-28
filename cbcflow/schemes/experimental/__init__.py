
### Experimental schemes (i.e. work in progress, not tested well, not documented)
# The simplest ipcs schemes with different optimizations
from .experimental.ipcs_opt_seg import SegregatedIPCS_Optimized

# Schemes with support for penalty pressure BCs
from .experimental.ipcs_penalty import PenaltyIPCS
from .experimental.ipcs_penalty_segregated import SegregatedPenaltyIPCS

# DG schemes not working in parallell
from .experimental.bottipietro import BottiPietro
from .experimental.piso import PISO
from .experimental.karper import Karper

# Coupled schemes
from .experimental.couplednonlinear import CoupledNonLinear
from .experimental.coupled_picard import CoupledPicard
from .experimental.stokes import Stokes
from .experimental.coupledpreconditioned import CoupledPreconditoned
from .experimental.coupledpreconditioned_kam import CoupledPreconditonedKAM

# Collect all schemes in list automatically
from ...core.nsscheme import NSScheme
experimental_schemes = [v for v in globals().values()
                        if hasattr(v, 'mro')
                        and issubclass(v, NSScheme)
                        and v is not NSScheme]
