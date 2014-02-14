from strip_code import strip_code
from spaces import NSSpacePool, NSSpacePoolMixed, NSSpacePoolSegregated, NSSpacePoolSplit, SpacePool
from show import show_problem

__all__ = [k for k,v in globals().items()
           if hasattr(v, "__module__")
           and __package__ in v.__module__]