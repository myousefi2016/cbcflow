#!/usr/bin/env python
__author__ = "Martin Sandve Alnaes <martinal@simula.no>"
__date__ = "2013-04-24"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from headflow import *
from headflow.dol import *

from nsdaproblem import NSDAProblem

class Problem(NSDAProblem):
    "3D pipe optimization test problem with known analytical solution."

    def __init__(self, params=None):
        NSDAProblem.__init__(self, params)

        # TODO: Need to fix parameter system to add parameters in an additional subclass level
        self.params.meshres = 32 # 64
        res = self.params.meshres

        # ... Problem specific definitions assumed to be present:

        mesh_filename = "mesh_r%d.xml.gz" % res
        mesh = Mesh(mesh_filename)

        facet_domains_filename = "facet_domains_r%d.xml.gz" % res
        facet_domains = MeshFunction("size_t", mesh, facet_domains_filename)

        self.initialize_geometry(mesh, facet_domains=facet_domains)

        self.wall_boundaries = (0,1)
        self.given_pressure_boundaries = (2,)
        self.control_boundaries = (3,4)

    def observations(self, spaces, t):
        "Return a list of observation functions that may need updating each timestep."
        U = spaces.U
        d = spaces.d

        # FIXME: Read a list of (tk,zk) tuples from file, zk = z at tk
        z = []
        num_z = 5
        ze = Expression("0.3 + 0.7*pow(sin(2*pi*t/period),2)", pi=DOLFIN_PI, t=0.0, period=1.0)
        ze.period = self.params.period
        for k in range(num_z):
            tk = self.params.T0 + (self.params.T-self.params.T0)*(k+1.0)/(num_z+2.0)

            # Functions to hold the value of z at the current timestep
            zk = as_vector([Function(U, name="z_%d_%d" % (k,i)) for i in xrange(d)])
            ze.t = tk
            zk[1].interpolate(ze)

            z.append((tk, zk))

        # Return observations tuple
        observations = (z,)
        return observations

if __name__ == "__main__":
    p = Problem()
    show_problem(p)
