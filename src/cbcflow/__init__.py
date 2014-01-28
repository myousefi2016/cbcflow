__version__ = "0.1"
__date__ = "2013-11-26"
__license__  = "GNU GPL version 3 or any later version"

# Basic utilities
from .core.paramdict import ParamDict

# Core component interfaces
from .core.nsproblem import NSProblem
from .core.nspostprocessor import NSPostProcessor
from .core.nsscheme import NSScheme
from .core.nssolver import NSSolver
from .core.nsreplay import NSReplay

# Problem inspection utilities
from .core.show import show_problem

# Timestepping utilities
from .core.adaptivetimestepping import AdaptiveTimestepping
from .core.constanttimestepping import ConstantTimestepping

# Boundary condition utilities
from .bcs.Pouseille import Pouseille, make_pouseille_bcs
from .bcs.Womersley import Womersley, make_womersley_bcs
from .bcs.Resistance import Resistance
from .bcs.UniformShear import UniformShear

# Postprocessing utilities
from .fields import show_fields, all_fields
for f in all_fields:
    exec("from .fields import %s" % (f,))

# Import schemes
from .schemes import IPCS
from .schemes import SegregatedIPCS
from .schemes import SegregatedIPCS_Optimized
from .schemes import IPCS_Stable
from .schemes import IPCS_Stabilized
from .schemes import PISO
from .schemes import Karper
from .schemes import BottiPietro
from .schemes import CoupledNonLinear
from .schemes import PenaltyIPCS
from .schemes import SegregatedPenaltyIPCS
from .schemes import CoupledPicard
from .schemes import Stokes
from .schemes import CoupledPreconditoned
from .schemes import Yosida

from .schemes import all_schemes

# Import utils
from .utils.fenicstools import *
#from .utils.create_slicemesh import create_slicemesh
from .utils.Slice import Slice
