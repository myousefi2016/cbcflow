from __future__ import division
__author__ = "Martin Alnaes <martinal@simula.no> and Oeyvind Evju <oyvinev@simula.no>"
__date__ = "2013-04-26"
__copyright__ = "Copyright (C) 2013-2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from ..dol import *

from .paramdict import ParamDict
from .parameterized import Parameterized


class NSScheme(Parameterized):
    """Base class for all Navier-Stokes schemes.

    TODO: Clean up and document new interface.
    """

    def __init__(self, params=None):
        Parameterized.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = ParamDict(
            # Discretization parameters
            u_degree = None, # Require this to be set explicitly by scheme subclass!
            p_degree = None, # Require this to be set explicitly by scheme subclass!

            # TODO: Move to nssolver?
            enable_timer=False,

            # TODO: These ipcs solvers are scheme specific and should maybe not be here,
            #       however they are used by most of the splitting schemes...
            # TODO: Split these into separate parameters for solver/preconditioner?
            solver_u_tent=("gmres", "hypre_euclid"),
            solver_p_neumann=("gmres", "hypre_amg"),
            solver_p_dirichlet=("gmres", "ml_amg"),
            solver_p=None, # overrides neumann/dirichlet if given
            solver_u_corr=("bicgstab", "hypre_euclid"),

            # Timestepping method
            adaptive_timestepping=False,
            )
        return params

    def solve(self, problem, update):
        "Solve Navier-Stokes problem by executing scheme."
        raise NotImplementedError("Scheme must implement solve method!")