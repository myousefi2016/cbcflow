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
from .bcs.Poiseuille import Poiseuille, make_poiseuille_bcs
from .bcs.Womersley import Womersley, make_womersley_bcs
from .bcs.Resistance import Resistance
from .bcs.UniformShear import UniformShear

# Postprocessing utilities and fields
from .fields import *

# Navier-Stokes solver schemes
from .schemes import *

# Import utils
from .utils.fenicstools import *
#from .utils.create_slicemesh import create_slicemesh
from .utils.Slice import Slice
