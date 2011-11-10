from __future__ import division
from problembase import *
from scipy import *
from numpy import array
from math import pi
from scipy.interpolate import splrep, splev

class Inflow(Expression):

    def __init__(self, V, problem):
        self.V = V
        self.mesh = V.mesh()
        self.problem = problem
        self.t_period = 1

	t  = array([0.  ,    27.,    42.,    58.,    69.,    88.,   110.,   130.,
                    136.,   168.,   201.,   254.,   274.,   290.,   312.,   325.,
                    347.,   365.,   402.,   425.,   440.,   491.,   546.,   618.,
                    703.,   758.,   828.,   897.,  1002.]) / (75/60.0) / 1000

	scale = 750
	#Create interpolated mean velocity in time
	v = array([390.        ,  398.76132931,  512.65861027,  642.32628399,
                   710.66465257,  770.24169184,  779.00302115,  817.55287009,
                   877.12990937,  941.96374622,  970.        ,  961.2386707 ,
                   910.42296073,  870.12084592,  843.83685801,  794.7734139 ,
                   694.89425982,  714.16918429,  682.62839879,  644.07854985,
                   647.58308157,  589.75830816,  559.96978852,  516.16314199,
                   486.37462236,  474.10876133,  456.58610272,  432.05438066,
                   390.]) / 574.211239628 * scale

        self.inflow = splrep(t, v)

    def eval_cell(self, values, x, ufc_cell):
        cell = Cell(self.mesh, ufc_cell.index)
        n = cell.normal(ufc_cell.local_facet)

	t = self.problem.t
	val = splev(t - int(t/self.t_period)*self.t_period, self.inflow)
	values[0] = -n.x()*val
	values[1] = -n.y()*val
	values[2] = -n.z()*val

    def value_shape(self):
        return (3,)

class InflowComp(Expression):

    def __init__(self, V, problem, comp):
        self.V = V
        self.mesh = V.mesh()
        self.problem = problem
        self.t_period = 1
        self.comp = comp

	t  = array([0.  ,    27.,    42.,    58.,    69.,    88.,   110.,   130.,
                    136.,   168.,   201.,   254.,   274.,   290.,   312.,   325.,
                    347.,   365.,   402.,   425.,   440.,   491.,   546.,   618.,
                    703.,   758.,   828.,   897.,  1002.]) / (75/60.0) / 1000

	scale = 750
	#Create interpolated mean velocity in time
	v = array([390.        ,  398.76132931,  512.65861027,  642.32628399,
                   710.66465257,  770.24169184,  779.00302115,  817.55287009,
                   877.12990937,  941.96374622,  970.        ,  961.2386707 ,
                   910.42296073,  870.12084592,  843.83685801,  794.7734139 ,
                   694.89425982,  714.16918429,  682.62839879,  644.07854985,
                   647.58308157,  589.75830816,  559.96978852,  516.16314199,
                   486.37462236,  474.10876133,  456.58610272,  432.05438066,
                   390.]) / 574.211239628 * scale


        self.inflow = splrep(t, v)

    def eval_cell(self, values, x, ufc_cell):
        cell = Cell(self.mesh, ufc_cell.index)
        n = cell.normal(ufc_cell.local_facet)

	t = self.problem.t
	val = splev(t - int(t/self.t_period)*self.t_period, self.inflow)
        comp = self.comp
        if comp == 0:
	    values[0] = -n.x()*val
        elif comp == 1:
	    values[0] = -n.y()*val
        elif comp == 2:
	    values[0] = -n.z()*val


# Define symmetric gradient
def epsilon(u):
    return grad(u)
#    return 0.5*(grad(u) + (grad(u).T))

# Problem definition
class Problem(ProblemBase):
    "3D artery with a saccular aneurysm."

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Load mesh
        refinement_level = options["refinement_level"]
        if refinement_level > 4:
            raise RuntimeError, "No mesh available for refinement level %d" % refinement_level
        self.mesh = Mesh("data/Aneurysm.xml.gz")

        # The body force term
        self.f = Constant((0, 0, 0))

        # Set viscosity
        self.nu = 3.5
        self.U = 1000

        # Set end-time
        self.T = 0.05
        self.First = True

    def initial_conditions(self, V, Q):
        if self.options['segregated']:
            zero = Constant(0)
            return (zero, zero, zero), zero
        else:
            return Constant((0, 0, 0)), Constant(0)

    def boundary_conditions(self, V, Q, t):
        if self.options['segregated']:
	    # Create no-slip boundary condition for velocity
	    self.g_noslip = Constant(0)
            bc_noslip = [DirichletBC(V, self.g_noslip, 0) for d in range(3)]

	     # Create inflow boundary condition for velocity
            self.g_inflow = [InflowComp(V, self, d) for d in range(3)]
            bc_inflow = [DirichletBC(V, g, 1) for g in self.g_inflow]

            bc_u = zip(bc_noslip, bc_inflow)

        else:
	    # Create no-slip boundary condition for velocity
	    self.g_noslip = Constant((0, 0, 0))
	    bc_noslip = DirichletBC(V, self.g_noslip, 0)

	     # Create inflow boundary condition for velocity
	    self.g_inflow = Inflow(V, self)
	    bc_inflow = DirichletBC(V, self.g_inflow, 1)

            bc_u = (bc_noslip, bc_inflow)

        # Create outflow boundary condition for pressure
        self.g_outflow = Constant(0)
        bc_outflow = [DirichletBC(Q, self.g_outflow, subd) for subd in (2,3)]
        bc_p = bc_outflow

        return bc_u, bc_p

    def pressure_bc(self, Q):
        return 0

    def update(self, t, u, p):
        self.t = t

    def functional(self, t, u, p):
         if t < self.T:
             return 0.0

         x = (0.025, -0.006, 0.0)
         if self.options['segregated']:
             return u[0](x)
         else:
             return u(x)[0]

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
