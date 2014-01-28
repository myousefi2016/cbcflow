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
    """Base class for all Navier-Stokes problems."""

    def __init__(self, params):
        """Initialize problem instance.

        Params will be taken from default_params and overridden
        by the values given to this constructor.
        """
        Parameterized.__init__(self, params)

        # Optionally set end time based on period and number of periods
        if self.params.T is None:
            if (self.params.num_periods is None or self.params.period is None):
                raise ValueError("You must provide parameter values for either end time T, or period and num_periods.")
            self.params.T = self.T0 + self.params.period * self.params.num_periods
        else:
            if self.params.num_periods is not None:
                raise ValueError("Ambiguous time period, cannot set both T and num_periods.")

    @classmethod
    def default_params(cls):
        """Returns the default parameters for a problem.

        Explanation of parameters:

        Time parameters:

          - start_timestep: int, initial time step number
          - dt: float, time discretization value
          - T0: float, initial time
          - T: float, end time
          - period: float, length of period
          - num_periods: float, number of periods to run

        Either T or period and num_period must be set.
        If T is not set, T=T0+period*num_periods is used.

        Physical parameters:

          - mu: float, kinematic viscosity
          - rho: float, mass density
            
        Space discretization parameters:

          - mesh_file: str, filename to load mesh from (if any)

        """
        params = ParamDict(
            # Physical parameters:
            mu=None,
            rho=None,

            # Time parameters:
            start_timestep = 0,
            dt=None,
            T0=0.0,
            T=None,
            period=None, 
            num_periods=None,

            # Spatial discretization parameters:
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

        Optimization problem support is currently experimental.
        Can be ignored for non-control problems.
        TODO: Document expected observations behaviour here.
        """
        return []

    def controls(self, spaces):
        """Return controls for optimization problem.

        Optimization problem support is currently experimental.
        Can be ignored for non-control problems.
        TODO: Document expected controls behaviour here.
        """
        return []

    def initial_conditions(self, spaces, controls):
        """Return initial conditions.

        TODO: Document expected initial_conditions behaviour here.

        Returns: u, p
        """
        raise NotImplementedError("initial_conditions must be overridden in subclass")

    def boundary_conditions(self, spaces, u, p, t, controls):
        """Return boundary conditions in raw format.

        TODO: Document expected boundary_conditions behaviour here.
        """
        raise NotImplementedError("boundary_conditions must be overridden in subclass")

    def body_force(self, spaces, t):
        """"Return body force, defaults to 0.

        TODO: Document expected body_force behaviour here.
        """
        d = self.mesh.geometry().dim()
        c0 = Constant(0.0)
        return [c0]*d

    def update(self, spaces, u, p, t, timestep, boundary_conditions, observations, controls): # TODO: Add body_force here
        """Update functions previously returned to new timestep.

        TODO: Document expected update behaviour here.

        The arguments bcs, observations, controls should be the
        exact lists of objects returned by boundary_conditions, observations, controls
        """
        pass

    def analytical_solution(self, spaces, t):
        """Return analytical solution.

        Can be ignored when no such solution exists,
        this is only used in the validation frameworks to
        validate schemes and test grid convergence etc.

        TODO: Document expected analytical_solution behaviour here.

        Returns: u, p
        """
        raise NotImplementedError("analytical_solution must be overridden in problem subclass to use analytical solution fields")

    def test_functionals(self, spaces):
        """Return fields to be used by regression tests.

        Can be ignored when no such solution exists,
        this is only used in the validation frameworks to
        validate schemes and test grid convergence etc.

        Returns: list of fields.
        """
        return []

    def test_references(self):
        """Return reference values corresponding to test_functionals to be used by regression tests.

        Can be ignored when no such solution exists,
        this is only used in the validation frameworks to
        validate schemes and test grid convergence etc.

        Returns: list of reference values.
        """
        return []
    
