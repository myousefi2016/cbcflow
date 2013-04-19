from __future__ import division

__author__ = "Joachim B Haga <jobh@simula.no>"
__date__ = "2011-11-11"
__copyright__ = "Copyright (C) 2011 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Martin Alnaes, 2013.

from headflow.problembase import *
from scipy import *
from numpy import array
from math import pi
from scipy.interpolate import splrep, splev

class InflowData(object):
    scale = 750

    def __init__(self, V, problem):
        self.mesh = V.mesh()
        self.problem = problem
        self.t_period = 1

        t  = array([0.  ,    27.,    42.,    58.,    69.,    88.,   110.,   130.,
                    136.,   168.,   201.,   254.,   274.,   290.,   312.,   325.,
                    347.,   365.,   402.,   425.,   440.,   491.,   546.,   618.,
                    703.,   758.,   828.,   897.,  1002.]) / (75/60.0) / 1000

        # Create interpolated mean velocity in time
        v = array([390.        ,  398.76132931,  512.65861027,  642.32628399,
                   710.66465257,  770.24169184,  779.00302115,  817.55287009,
                   877.12990937,  941.96374622,  970.        ,  961.2386707 ,
                   910.42296073,  870.12084592,  843.83685801,  794.7734139 ,
                   694.89425982,  714.16918429,  682.62839879,  644.07854985,
                   647.58308157,  589.75830816,  559.96978852,  516.16314199,
                   486.37462236,  474.10876133,  456.58610272,  432.05438066,
                   390.]) / 574.211239628 * self.scale

        self.inflow_spline = splrep(t, v)

    def __call__(self, x, ufc_cell):
        val = splev(self.problem.t % self.t_period, self.inflow_spline)

        cell = Cell(self.mesh, ufc_cell.index)
        n = cell.normal(ufc_cell.local_facet)
        return [-n.x()*val, -n.y()*val, -n.z()*val]

class InflowVec(Expression):
    def __init__(self, V, problem):
        self.data = InflowData(V, problem)
    def eval_cell(self, values, x, ufc_cell):
        values[:] = self.data(x, ufc_cell)
    def value_shape(self):
        return 3,

class InflowComp(Expression):
    def __init__(self, V, problem, component):
        self.data = InflowData(V, problem)
        self.component = component
    def eval_cell(self, values, x, ufc_cell):
        values[0] = self.data(x, ufc_cell)[self.component]

# Problem definition
class Problem(NSProblem):
    "3D artery with a saccular aneurysm."

    def __init__(self, params):
        NSProblem.__init__(self, params)

        # Load mesh
        self.mesh = Mesh(retrieve("data/Aneurysm.xml.gz"))

        # The body force term
        self.f = self.uConstant((0, 0, 0))

        # Set viscosity
        self.nu = 3.45

        # Characteristic velocity in the domain (used to determine timestep)
        self.U = InflowData.scale

        # Set current and end-time
        self.t = 0.0
        self.T = 0.5

    def initial_conditions(self, V, Q):
        return self.uConstant((0, 0, 0)) + [Constant(0)]

    def resistance(self, mesh, mark, C, p0):
        if not hasattr(self, "u"):
            if master:
                print "self.u not initialized, assuming zero flux (resistance is %.3g)"%p0
            return Constant(p0)
        n = FacetNormal(mesh)
        flux = inner(self.u, n)*ds(mark)
        Q = assemble(flux, mesh=mesh)
        R = C*Q + p0
        if master:
            print "Computed resistance over marker %d is %.3g, the flux is %.3g"%(mark, R, Q)
        return Constant(R)

    def boundary_conditions(self, V, Q, t):
        # Create no-slip boundary condition for velocity
        self.g_noslip = self.uConstant((0, 0, 0))
        bc_noslip = [DirichletBC(V, g, 0) for g in self.g_noslip]

        # Create inflow boundary condition for velocity
        if self.params.segregated:
            self.g_inflow = [InflowComp(V, self, d) for d in range(3)]
        else:
            self.g_inflow = [InflowVec(V, self)]
        bc_inflow = [DirichletBC(V, g, 1) for g in self.g_inflow]

        # Create outflow boundary condition for pressure
        self.g_outflow = {}
        self.g_outflow[2] = self.resistance(V.mesh(), 2, C=5.97, p0=0)
        self.g_outflow[3] = self.resistance(V.mesh(), 3, C=5.97, p0=0)
        bc_outflow = [DirichletBC(Q, self.g_outflow[marker], marker) for marker in (2,3)]

        bc_u = zip(bc_inflow, bc_noslip) # Important: inflow before noslip
        bc_p = [bc_outflow]

        return bc_u + bc_p

    def update(self, t, u, p):
        self.t = t

    def functional(self, t, u, p):
         if t < self.T:
             return 0.0

         return self.uEval(u, 0, (0.025, -0.006, 0.0))

    def reference(self, t):
        """The reference value was computed using on a fine mesh
        (level 6). Values obtained for different refinement levels
        are listed below for Chorin and IPCS.

              Chorin                 IPCS
        ----------------------------------------
        -0.0325040608617000  -0.0333250879034000
        -0.0470001557641000  -0.0458749339862000
        -0.0370348732066000  -0.0364138324117000
        -0.0359768558469000  -0.0358236703894000
        -0.0356064894317000  -0.0354277722246000
        -0.0355250220872000  -0.0353312047875000
        -0.0356105862451000  -0.0354251625379000

        The reference value is taken as the average of the values
        for Chorin and IPCS on the finest mesh.
        """
        if t < self.T:
            return 0.0

        return -0.0355

    def __str__(self):
        return "Aneurysm"

