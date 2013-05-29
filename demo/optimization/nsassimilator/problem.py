#!/usr/bin/env python
__author__ = "Martin Sandve Alnaes <martinal@simula.no>"
__date__ = "2013-04-24"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from headflow import *
from headflow.dol import *

from daproblem import NSDAProblem

class Problem(NSDAProblem):
    "3D pipe optimization test problem with known analytical solution."

    def __init__(self, params=None):
        NSDAProblem.__init__(self, params)

        # ... Problem specific definitions assumed to be present:

        # Get 3D pipe mesh from file
        mesh = Mesh("../../../data/pipe_0.2.xml.gz")
        self.initialize_geometry(mesh)

        # Boundary markers for different parts
        self.wall_boundaries = (0,)
        self.given_pressure_boundaries = (2,)
        self.control_boundaries = (1,)

        # ... Purely problem specific variables:

        # Known properties of the mesh
        self.length = 10.0
        self.radius = 0.5
        self.beta = 4*(self.params.mu/self.params.rho) / self.radius**2

    def x_observation_expression(self, t):
        # Quadratic profile times a transient pulse
        ze = "(minflow + (1.0-minflow)*pow(sin(2.0*DOLFIN_PI*t/period),2)) * upeak * (r*r-x[1]*x[1]-x[2]*x[2]) / (r*r)"
        ze = Expression(ze, upeak=1.0, r=0.5, period=1.0, minflow=0.3, t=0.0, name="ze")
        ze.r = self.radius
        ze.period = self.params.period
        ze.minflow = 0.3
        ze.upeak = 1.0
        ze.t = float(t)
        return ze

    def _observation_expression(self, t):
        ze = "upeak * (r*r-x[1]*x[1]-x[2]*x[2]) / (r*r)"
        ze = Expression(ze, upeak=1.0, r=0.5, name="ze")
        return ze

    def observations(self, spaces, t):
        "Return a list of observation functions that may need updating each timestep."
        U = spaces.U
        d = spaces.d


        # Function to hold the value of z at the current timestep
        z = as_vector([Function(U, name="z_%d" % i) for i in xrange(d)])
        ze = self._observation_expression(t)
        z[0].interpolate(ze)

        # FIXME: Return a list of (tk,zk) tuples, zk = z at tk instead


        # Return observations tuple
        observations = (z,)
        return observations

    def update_observations(self, spaces, t, observations):
        "Update functions in list returned by auxilliary_functions() for this timestep."
        z, = observations

        # Interpolate ze into x-component of z
        ze = self._observation_expression(t)
        z[0].interpolate(ze)

    def initial_control_values(self):
        d = self.mesh.ufl_cell().d
        e0 = Expression("0.0")
        ze = "upeak * (r*r-x[1]*x[1]-x[2]*x[2]) / (r*r)"
        ze = Expression(ze, upeak=1.0, r=0.5, name="ze")
        u0 = [ze, e0, e0][:d]

        # Start with a fraction of the actual solution pressure drop to get some movement
        num_p_controls = len(self.control_boundaries)
        p_coeffs0 = [-0.1*self.length*self.beta]*num_p_controls

        return u0, p_coeffs0


if __name__ == "__main__":
    p = Problem()
    show_problem(p)
