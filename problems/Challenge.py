from __future__ import division
from problembase import *
from scipy import *
from numpy import array, sin, cos, pi
from scipy.interpolate import splrep, splev

"""
SI units:
N = kg m / s2
Pa = N / m2 = kg / (m s2)
1 ml = 1 cm3

cgs system:
1 dyn = 1 g cm / s2 = 10−5 N
1 dyn / cm2 = 1 g / (cm s2) = 1e-3 kg / (1e-2 m s2) = 1e-1 Pa
1 poise = 1 P = 1 g / (cm s) = 10−1 Pa s
1 cP = 0.01 P = 1 mPa s = 0.001 Pa s

Units from challenge readme:
mu = 0.04 Poise = 0.04 g / (cm s)
rho = 1.0 g / cm3
Q = ml / s = cm3 / s
nu = mu / rho = 0.04 (g / (cm s)) / (g / cm3) = 0.04 (cm2 / s)

Checking units:
[Dv/dt] = [nu lap u] = [p/rho]
[Dv/dt] = cm/s2
[nu laplace u] = cm2/s  1/cm2  cm/s = cm/s2
[grad p/rho] = [p] cm2/g = cm/s2,  [p] = g/(cm s2) = dyn

Area of inlet is approx. 0.25 cm2, computed from A = pi * (0.56 cm / 2)**2
Average inflow velocity is U = Q / A ~ (5 cm3/s) / (0.25 cm2) = 20 cm/s

EXPERIMENTS
===========

PULSATILE
---------
Two different pulsatile flows, both based on the same waveform
shape, but with different mean flow rates:
Qpulsatile = (5.13, 6.41) # cm3 / s

These two pulsatile cases correspond to mean WSS at inlet of
wss_pulsatile = (12.0, 15.0) # 1 dyne/cm2 = 0.1 Pa

STATIONARY
----------
Four steady flows at rates:
Qsteady = (5.13, 6.41, 9.14, 11.42) # cm3 / s

Note that the first two test cases are the same as the pulsatile flow rates.

"""


A = [1, -0.23313344, -0.11235758, 0.10141715, 0.06681337, -0.044572343, -0.055327477, 0.040199067, 0.01279207, -0.002555173, -0.006805238, 0.002761498, -0.003147682, 0.003569664, 0.005402948, -0.002816467, 0.000163798, 8.38311E-05, -0.001517142, 0.001394522, 0.00044339, -0.000565792, -6.48123E-05] 

B = [0, 0.145238823, -0.095805132, -0.117147521, 0.07563348, 0.060636658, -0.046028338, -0.031658495, 0.015095811, 0.01114202, 0.001937877, -0.003619434, 0.000382924, -0.005482582, 0.003510867, 0.003397822, -0.000521362, 0.000866551, -0.001248326, -0.00076668, 0.001208502, 0.000163361, 0.000388013]


counter = 0 
N = 100 

def time_dependent_velocity(t):
    "Returns unitless temporal variation of velocity."
    velocity = 0 
    for k in range(len(A)): 
        velocity += A[k]*cos(2*pi*k*t)
        velocity += B[k]*sin(2*pi*k*t)
    return velocity

class InflowData(object):

    def __init__(self, V, problem):
        self.mesh = V.mesh()
        self.problem = problem
        self.velocity = problem.velocity
        self.val = self.velocity 
        self.stationary = problem.stationary
        self.t = 0 
        self.N = 100 

    def __call__(self, x, ufc_cell):
        cell = Cell(self.mesh, ufc_cell.index)
        n = cell.normal(ufc_cell.local_facet)

        global counter
        global N

        if self.problem.t > self.t and not self.problem.stationary:
            self.t = self.problem.t
            self.val = self.velocity*time_dependent_velocity(self.t)
        if self.problem.stationary and counter <= N:
            if self.problem.t > self.t:
                self.t = self.problem.t
                counter += 1
                self.val = float(self.velocity*self.counter)/self.N
                print self.val, self.velocity, self.counter, self.N

        val = self.val
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
        refinement = self.options["refinement_level"]
        self.stationary = self.options["stationary"]
        boundary_layers = self.options["boundary_layers"]

        if boundary_layers:
            mesh_filename = {
                0: "mesh_750k_BL_t.xml.gz",
                1: "mesh_2mio_BL_t.xml.gz",
                2: "mesh_4mio_BL_t.xml.gz",
                }[refinement]
        else:
            mesh_filename = {
                0: "mesh_500k.xml.gz",
                1: "mesh_1mio.xml.gz",
                2: "mesh_2mio.xml.gz",
                3: "mesh_4mio.xml.gz",
            }[refinement]
        self.mesh = Mesh(os.path.join("data", "challenge", mesh_filename))

        # Select flux in cm3/s
        self.testcase = self.options["test_case"]
        self.flux = {
            1: 5.13,
            2: 6.41,
            3: 9.14,
            4: 11.42,
            }[self.testcase]

        # The body force term
        self.f = self.uConstant((0, 0, 0))

        # Set viscosity (cm2/s)
        self.nu = 0.04

        # Set current and end-time
        self.t = 0.0
        self.T = 0.05 # s # TODO: Is this an input parameter somewhere? Add parameter if not! Peak systole is 0.275 s.

        # Compute volume (cm3) and areas (cm2)
        one = Constant(1)
        self.V0 = assemble(one*dx, mesh=self.mesh)
        self.A0 = assemble(one*ds(0), mesh=self.mesh)
        self.A1 = assemble(one*ds(1), mesh=self.mesh)
        self.A2 = assemble(one*ds(2), mesh=self.mesh)

        print "Volume of the geometry is (dx)   ", self.V0, "cm^3"
        print "Areal  of the no-slip is (ds(0)  ", self.A0, "cm^2"
        print "Areal  of the inflow is (ds(1))  ", self.A1, "cm^2"
        print "Areal  of the outflow is (ds(2)) ", self.A2, "cm^2"

        # cm/s
        self.velocity = self.flux / self.A1 

        # Characteristic velocity (U) in the domain (used to determine timestep)
        cfl_factor = 16 # TODO: Is this a bit high? Reduce by a factor 2-3 to speed up?
        self.U = self.velocity*cfl_factor
        h = MPI.min(self.mesh.hmin())
        print "Characteristic velocity set to", self.U, "cm/s"
        print "mesh size          ", h, "cm"
        print "velocity at inflow ", self.velocity, "cm/s"
        print "Number of cells    ", self.mesh.num_cells()
        print "Number of vertices ", self.mesh.num_vertices()

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
        bc_outflow = [DirichletBC(Q, self.g_outflow, marker) for marker in (2,)]

        bc_u = zip(bc_inflow, bc_noslip) # Important: inflow before noslip
        bc_p = [bc_outflow]

        return bc_u + bc_p

    def pressure_bc(self, Q):
        return 0

    def update(self, t, u, p):
        self.t = t

    def functional(self, t, u, p):
        n = FacetNormal(self.mesh)
        b0 = assemble(dot(u,n)*ds(0)) # cm2/s
        b1 = assemble(dot(u,n)*ds(1)) 
        b2 = assemble(dot(u,n)*ds(2)) 
        b3 = assemble(dot(u,n)*ds(3))

        #FIXME should use selected points
        p_max = p.vector().max()
        p_min = p.vector().min()

        print "flux ds0 ", b0, "cm^2/s"
        print "flux ds1 ", b1, "cm^2/s"
        print "flux ds2 ", b2, "cm^2/s"
        print "flux ds3 ", b3, "cm^2/s"
        print "p_min ", p_min, "g/(cm s^2) or dyne"
        print "p_max ", p_max, "g/(cm s^2) or dyne"
        if self.options["segregated"]: 
            u_max = max(ui.vector().norm('linf') for ui in u) 
        else:
            u_max = u.vector().norm('linf')  
        print "u_max ", u_max, " U ", self.U, "cm/s"

        return p_max - p_min 

    def reference(self, t):
        return 0 

    def __str__(self):
        return "Challenge"

