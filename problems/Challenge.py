from __future__ import division
from problembase import *
from scipy import *
from numpy import array
from math import pi
from scipy.interpolate import splrep, splev

class InflowData(object):
    scale = 40 

    def __init__(self, V, problem):
        self.mesh = V.mesh()
        self.problem = problem
        self.val = problem.velocity 


    def __call__(self, x, ufc_cell):
        val = self.val 
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
class Problem(ProblemBase):
    "3D artery with a saccular aneurysm."

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Load mesh
        self.mesh = Mesh("data/mesh_500k.xml.gz")

        self.refinement = self.options["refinement_level"] 
        self.testcase = self.options["test_case"] 
        self.flux = 0
        if self.testcase == 1: 
	    self.flux = 5.13 
        elif self.testcase == 2:
	    self.flux = 6.41 
        elif self.testcase == 3:
	    self.flux = 9.14 
        elif self.testcase == 4:
	    self.flux = 11.42 


        # The body force term
        self.f = self.uConstant((0, 0, 0))

        # Set viscosity
        self.nu = 0.04 


        # Set current and end-time
        self.t = 0.0
        self.T = 0.5


        one = Constant(1)
	self.V0 = assemble(one*dx, mesh=self.mesh)
	self.A0 = assemble(one*ds(0), mesh=self.mesh)
	self.A1 = assemble(one*ds(1), mesh=self.mesh)
	self.A2 = assemble(one*ds(2), mesh=self.mesh)

	print "Volume of the geometry is ", self.V0 
	print "Areal  of the no-slip is  ", self.A0 
	print "Areal  of the inflow is   ", self.A1 
	print "Areal  of the outflow is  ", self.A2 

	self.velocity = self.flux / self.A1 

        # Characteristic velocity in the domain (used to determine timestep)
        self.U = self.velocity*4  


    def initial_conditions(self, V, Q):
        return self.uConstant((0, 0, 0)) + [Constant(0)]

    def boundary_conditions(self, V, Q, t):
        # Create no-slip boundary condition for velocity
        self.g_noslip = self.uConstant((0, 0, 0))
        bc_noslip = [DirichletBC(V, g, 0) for g in self.g_noslip]

        # Create inflow boundary condition for velocity
        if self.options['segregated']:
            self.g_inflow = [InflowComp(V, self, d) for d in range(3)]
        else:
            self.g_inflow = [InflowVec(V, self)]
        bc_inflow = [DirichletBC(V, g, 1) for g in self.g_inflow]

        # Create outflow boundary condition for pressure
        self.g_outflow = Constant(0)
        bc_outflow = [DirichletBC(Q, self.g_outflow, marker) for marker in (2,3)]

        bc_u = zip(bc_inflow, bc_noslip) # Important: inflow before noslip
        bc_p = [bc_outflow]

        return bc_u + bc_p

    def pressure_bc(self, Q):
        return 0

    def update(self, t, u, p):
        self.t = t

    def functional(self, t, u, p):

         n = FacetNormal(self.mesh)
         b = dot(u,n)*ds(1) 
         p_max = p.vector().max()
         p_min = p.vector().min()

         print "flux ", b 
         print "p_min ", p_min 
         print "p_max ", p_max
         if self.options["segregated"]: 
             u_max = max(ui.vector().norm('linf') for ui in u) 
         else:
             u_max = u[0].vector().norm('linf')  

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
