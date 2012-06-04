from __future__ import division

__author__ = "Joachim B Haga <jobh@simula.no>"
__date__ = "2011-11-11"
__copyright__ = "Copyright (C) 2011 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from problembase import *
from scipy import *
from numpy import array
from math import pi
from scipy.interpolate import splrep, splev

class InflowData(object):

    def __init__(self, V, t, problem, inlet):
        self.mesh = V.mesh()
        self.center = problem.center[inlet]
        self.radius = problem.radius[inlet]
        self.peakU = problem.peakU[inlet]
        self.t = t
        
    #def __call__(self, x, ufc_cell):
    #    val = 1
    #    cell = Cell(self.mesh, ufc_cell.index)
    #    n = cell.normal(ufc_cell.local_facet)
    #    return [-n.x()*val, -n.y()*val, -n.z()*val]

    def __call__(self, x, ufc_cell):
        # print 'self.center is', self.center
	 r = sqrt( (self.center[0] - x[0])**2 + (self.center[1] - x[1])**2 + (self.center[2] - x[2])**2) 
         #pp = (1-(r/self.radius)**2)*(1-exp(-10*self.t))
         pp = (1-(r/self.radius)**2)
         #print "exp ", (1-exp(-10*self.t))         
         cell = Cell(self.mesh, ufc_cell.index) 
         n = cell.normal(ufc_cell.local_facet)
         return [-n.x()*self.peakU*pp, -n.y()*self.peakU*pp, -n.z()*self.peakU*pp]

class InflowVec(Expression):
    def __init__(self, V, problem, **kwargs):
        t = kwargs['t']
        inlet = kwargs['inlet'] 
        self.data = InflowData(V, t, problem, inlet)
    def eval_cell(self, values, x, ufc_cell):
        values[:] = self.data(x, ufc_cell)
    def value_shape(self):
        return 3,

class InflowComp(Expression):
    def __init__(self, V, problem, component, **kwargs):
        t = kwargs['t']
        inlet = kwargs['inlet']
        self.data = InflowData(V, t, problem, inlet)
        self.component = component
    def eval_cell(self, values, x, ufc_cell):
        values[0] = self.data(x, ufc_cell)[self.component]

# Problem definition
class Problem(ProblemBase):
    "3D T-junction with cylindrical branches."

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Load mesh
        self.mesh = Mesh(self.retrieve("data/Tjunction/DA90/mesh324k.xml"))
#        self.mesh = Mesh("data/Tjunction/DA90/mesh169k.xml.gz")
#        self.mesh = Mesh("data/Tjunction/DA90/mesh1066k.xml.gz")

        # The body force term
        self.f = self.uConstant((0, 0, 0))

        # Set viscosity
        self.nu = 0.035

        # Set peak velocity
        self.peakU = zeros(2)
        self.peakU[0] = 52
        self.peakU[1] = 78

        # Characteristic velocity in the domain (used to determine timestep)
        self.U = MPI.max(max(self.peakU))

        self.dimen = zeros(3)

        def modlen(dim):
            c = self.mesh.coordinates()
            return MPI.max(c[:,dim].max()) - MPI.min(c[:,dim].min())

        self.dimen[0] = modlen(0)
        print "Model length is ", self.dimen[0]

        self.dimen[1] = modlen(1)
        print "Model height is ", self.dimen[1]

        self.dimen[2] = modlen(2)
        print "Model depth is ", self.dimen[2]

        # Set current and end-time
        self.t = 0.0
        # calculate T based on mean (not peak) velocity
        self.T = self.options['T'] or 2.*self.dimen.max()/(mean(self.peakU)/2.)

        self.wall = 0
        self.inlet = zeros(2, 'int')
        self.inlet[0] = 1  # Left branch
        self.inlet[1] = 3  # Side branch
        self.outlet = 2  # Right branch
        self.radius = zeros(2)  # for each inlet
        self.center = zeros([2,3])
        self.A = zeros(2)

        self.pres = [0, 0]
        
        VV = VectorFunctionSpace(self.mesh, "CG",1)
        
        self.A[0] = assemble(Constant(1)*ds(self.inlet[0]), mesh=self.mesh)
        p = project(Expression(("x[0]", "x[1]", "x[2]")), VV)
        self.radius[0] = sqrt(self.A[0]/DOLFIN_PI)
        print 'Inflow1 radius is', self.radius[0]
        for i in range(3):
           self.center[0,i] = assemble(p[i]*ds(self.inlet[0]))/self.A[0]
        print 'Inflow1 center is ', self.center[0]

        self.A[1] = assemble(Constant(1)*ds(self.inlet[1]), mesh=self.mesh)
        p = project(Expression(("x[0]", "x[1]", "x[2]")), VV)
        self.radius[1] = sqrt(self.A[1]/DOLFIN_PI)
        print 'Inflow2 radius is', self.radius[1]
        for i in range(3):
            self.center[1,i] = assemble(p[i]*ds(self.inlet[1]))/self.A[1]
        print 'Inflow2 center is ', self.center[1]

    def initial_conditions(self, V, Q):
        return self.uConstant((0, 0, 0)) + [Constant(0)]

    def boundary_conditions(self, V, Q, t):

        # Create no-slip boundary condition for velocity
        self.g_noslip = self.uConstant((0, 0, 0))
        bc_noslip = [DirichletBC(V, g, self.wall) for g in self.g_noslip]

        # Create inflow boundary condition for velocity
        if self.options['segregated']:
            self.g_inflow1 = [InflowComp(V, self, d, t=t, inlet=0) for d in range(3)]
        else:
            self.g_inflow1 = [InflowVec(V, self, t=t, inlet=0)]        
        
        #self.g_inflow1 = self.uConstant((1, 0, 0))
        
        bc_inflow1 = [DirichletBC(V, g, self.inlet[0]) for g in self.g_inflow1]
	#bc_outflow1 = [DirichletBC(Q, 0, self.inlet1)]

        if self.options['segregated']:
            self.g_inflow2 = [InflowComp(V, self, d, t=t, inlet=1) for d in range(3)]
        else:
            self.g_inflow2 = [InflowVec(V, self, t=t, inlet=1)]

        # self.g_inflow2 = self.uConstant((1, 1, 1))
        bc_inflow2 = [DirichletBC(V, g, self.inlet[1]) for g in self.g_inflow2]
	#bc_outflow1 = [DirichletBC(Q, 0, self.inlet2)]
        
	# Create constant pressure boundary
	bc_outflow = [DirichletBC(Q, 0, self.outlet)]

        bc_u = zip(bc_inflow1, bc_inflow2, bc_noslip) # Important: inflow before noslip
        #bc_u = zip(bc_inflow2, bc_noslip) # Important: inflow before noslip
        #bc_p = zip(bc_outflow, bc_outflow1)
        bc_p = [bc_outflow]

        return bc_u + bc_p

    def update(self, t, u, p):
        self.t = t

    def functional(self, t, u, p):
         
        self.pp = p*ds(self.inlet[0])
        self.pres[0] = assemble(self.pp, mesh=self.mesh)/self.A[0]

        self.pp = p*ds(self.inlet[1])
        self.pres[1] = assemble(self.pp, mesh=self.mesh)/self.A[1]

        return self.pres

    def reference(self, t):
        return None

    def __str__(self):
        return "Tjunction"
