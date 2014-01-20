from __future__ import division
__author__ = "Martin Alnaes <martinal@simula.no> and Oeyvind Evju <oyvinev@simula.no>"
__date__ = "2013-04-26"
__copyright__ = "Copyright (C) 2013-2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

import ufl

from .paramdict import ParamDict
from .parameterized import Parameterized
from ..dol import Constant, MeshFunction

class NSProblem(Parameterized):
    """Base class for all Navier-Stokes problems.

    TODO: Document new interface.
    """

    def __init__(self, params):
        Parameterized.__init__(self, params)

        # Optionally set end time based on period and number of periods
        if self.params.T is None and (self.params.num_periods is None or self.params.period is None):
            raise RuntimeError("You must provide parameter values for either end time T, or period and num_periods.")
        if self.params.T is None:
            self.params.T = self.params.period * self.params.num_periods
        else:
            assert self.params.num_periods is None, "Cannot set both T and num_periods."

    @classmethod
    def default_params(cls):
        params = ParamDict(
            # Time parameters:
            start_timestep = 0,
            dt=None, # Timestep
            T0=0.0, # Initial time (default 0 is the usual choice)
            T=None, # End time
            period=None, 
            num_periods=None,

            # Physical parameters:
            mu=None,    # Kinematic viscosity
            rho=None,    # Density
            
            mesh_file=None,
            )
        return params

    def initialize_geometry(self, mesh, facet_domains=None, cell_domains=None):
        """Stores mesh, domains and related quantities in a canonical member naming.

        Creates properties:

            - mesh
            - facet_domains
            - cell_domains
            - ds
            - dS
            - dx

        """
        # Store geometry properties
        self.mesh = mesh
        self.facet_domains = facet_domains
        self.cell_domains = cell_domains

        # Fetch domains from mesh if necessary and avaiable
        domains = mesh.domains()
        if domains is not None:
            dim = mesh.geometry().dim()
            if self.facet_domains is None:
                self.facet_domains = MeshFunction("size_t", mesh, dim-1, domains)
            if self.cell_domains is None:
                self.cell_domains = MeshFunction("size_t", mesh, dim, domains)

        # Attach domains to measures for convenience
        self.ds = ufl.ds if self.facet_domains is None else ufl.ds[self.facet_domains]
        self.dS = ufl.dS if self.facet_domains is None else ufl.dS[self.facet_domains]
        self.dx = ufl.dx if self.cell_domains  is None else ufl.dx[self.cell_domains]

    def observations(self, spaces, t):
        """Return observations of velocity for optimization problem.

        Can be ignored for non-control problems.

        TODO: Document expected observations behaviour here.
        """
        return []

    def controls(self, spaces):
        """Return controls for optimization problem.

        Can be ignored for non-control problems.

        TODO: Document expected controls behaviour here.
        """
        return []

    def analytical_solution(self, spaces, t):
        """Return analytical solution.

        Can be ignored when no such solution exists,
        this is only used in the validation frameworks to
        validate schemes and test grid convergence etc.

        TODO: Document expected analytical_solution behaviour here.
        """
        raise NotImplementedError("analytical_solution must be overridden in subclass")

    def initial_conditions(self, spaces, controls):
        """Return initial conditions.

        TODO: Document expected initial_conditions behaviour here.
        """
        raise NotImplementedError("initial_conditions must be overridden in subclass")

    def boundary_conditions(self, spaces, u, p, t, controls):
        """Return boundary conditions in raw format.

        TODO: Document expected boundary_conditions behaviour here.
        """
        raise NotImplementedError("boundary_conditions must be overridden in subclass")

    def body_force(self, spaces, t):
        """"Return body force, defaults to 0.

        TODO: Document expected body_force behaviour here."""
        d = self.mesh.topology().dim()
        c0 = Constant(0.0)
        return [c0]*d

    def update(self, spaces, u, p, t, timestep, bcs, observations, controls):
        """Update functions previously returned to new timestep.

        TODO: Document expected update behaviour here.
        """
        pass
