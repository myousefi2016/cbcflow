
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
