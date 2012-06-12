
# Chiari test problem
#
# Kent-Andre Mardal 2008-05-27  
# changed Susanne Hentschel 2009-03-09
# Joachim Berdal Haga 2012

from problems.problembase import *
from numpy import array, pi

def get_pulse_input_function(V, z_index, factor, A, HR_inv, HR, b, f1):
    two_pi = 3.4 * pi
    rad = two_pi /HR_inv
    #v_z = "-factor*(-A*(exp(-fmod(t,T)*rad) * Ees * (sin( - f1*fmod(t,T)*rad) - vd) - (1-exp( - factor*fmod(t,T)*rad)) * p0 * (exp(sin( - fmod(t,T)*rad) - vd) -1) ) -b)"
    v_z = "-2.5 * factor"
    vel = ["0.0", "0.0", "0.0"]
    vel[z_index] = v_z

    return {'cppcode': vel, "factor":factor, "A":A, "p0":1, "vd":0.03, "Ees":50, "T": HR_inv, "HR":HR, "rad":rad, "b":b, 'f1':f1}

def get_sine_input_function(V, z_index, HR, HR_inv, v_max):
    v_z = "-sin(2*pi*HR*fmod(t,T))*(v_max)"
    #v_z = "-v_max"
    vel = ["0.0", "0.0", "0.0"]
    vel[z_index] = v_z
    return {'cppcode': vel, 'HR':HR, 'v_max':v_max, 'T':HR_inv}

# Problem definition
class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        case_name = options['casename'] or 'caseStenosis40'
        print 'use case ' + case_name
        exec("from cases import %s as case" %case_name)

        self.case = case()
        self.mesh = self.case.get_mesh()
        self.sub_domains = self.case.get_sub_domains()
        self.contour = self.case.get_contour()
        self.top = self.case.get_top()
        self.bottom = self.case.get_bottom()

        # define waveform
        self.waveform = options["waveform"]
        if self.waveform not in ["pulse", "sine"]:
            raise RuntimeError("Waveform '%s' not valid"%self.waveform)
        print "use waveform " + self.waveform

        # define pulse variable f1 for use in case of pulse
        if self.waveform == "pulse":
            self.f1 = options["f1"]
            print "use f1 = " + str(self.f1)

        # Set viscosity
        self.nu = 0.7 *10**(-2) #cm^2/s  ;    0.658*10**-6

        # Create right-hand side function
        self.f = self.uConstant((0.0, 0.0, 0.0))
        n = tetrahedron.n

        self.h = MPI.min(self.mesh.hmin())
        self.A0 = self.area(0)
        self.A1 = self.area(1)
        self.A2 = self.area(2)
        self.AS = 0.93  #reference
        self.D = sqrt(self.A1/pi) #characteristic diameter in cm
        self.z_index = self.case.get_z_index()

        #set heart rate
        self.HR = 1.16 #beats per second; from 70 bpm
        self.HR_inv = 1.0/self.HR

        # Set end-time
        self.T = options["T"] or 1.2

        self.dt = options["dt"] or 0.00025#0.01

        self.peak_v = 0.0
        self.pressure_at_peak_v = 0.0

#       self.pressure_index = 52000
#       self.pressure_index = 520
#       self.characteristic = "preassure at node "+ str(self.pressure_index) + " at peak velocity"
        

        # calculate maximum velocity and parameters to define input velocity function
        self.SV = 1.0#0.5 # gupta: 0.27 ml to and fro per cardiac cycle -Linge used 17 ml/s as peak -Martin 5.4 ml/s at peak -> 3.176 ml in (0.85/2) s -> SV = 3.176*sin_integration_factor*A1 = 0.976 ml (about 1 ml!!) up and down per cardiac cycle (A1 = 0.976) -However, the sin_integration_factor is for HR = 1 and the real SV used here is 0.85 ml!!!!!
        self.flow_per_unit_area1 = self.SV/self.A1
        self.flow_per_unit_area2 = self.SV/self.A2
    #   self.flow_per_unit_area = self.SV/self.AS
        # for pulse function:
        if self.waveform == "pulse":
            self.A = 2.9/16 # scale to get max = 2.5 and min = -0.5 for f1 = 1
            self.factor1 = self.flow_per_unit_area1/(0.324) #factor that scales the function so that systolic and diastolic flow equal volume flow  
            self.factor2 = self.flow_per_unit_area2/(0.324)
            self.v_max1 = 2.5 * self.factor1
            self.v_max2 = 2.5 * self.factor2
            self.b = 0.465  # translating the function "down"
        # for sine function
        elif self.waveform =="sine":
            sin_integration_factor = 0.315 #(integral of the positive part of a sine, scaled on an interval length of 1)
            self.v_max1 = self.flow_per_unit_area1/sin_integration_factor
            self.v_max2 = self.flow_per_unit_area2/sin_integration_factor
        #scale velocity so that stroke volume is equal to volume flow

        self.print_info()
        #self.get_characteristic_value()


    def get_characteristic(self):
        return self.characteristic

    def get_characteristic_value(self):
        return self.pressure_at_peak_v

    def print_info(self):
        print "h_min ", self.h
        print "dt ", self.dt
        print "A0 =", self.A0, "(no-slip)"
        print "A1 =", self.A1, "(top)"
        print "A2 =", self.A2, "(bottom)"
        print "suggested cross section A =", self.AS, "cm^2"
        print "v_max on top    ", self.v_max1
        print "v_max bottom    ", self.v_max2
        print "Re (top) = ", self.v_max1 * self.D / self.nu
        print "Re (bottom) = ", self.v_max2 * self.D / self.nu

    def update(self, t, u, p):
        for g in self.g0 + self.g1 + self.g2:
            g.t = t

    def boundary_conditions(self, V, Q, t):

        # Create no-slip boundary condition for velocity
        self.g0 = self.uConstant((0.0, 0.0, 0.0))
        #bc0 = DirichletBC(V, self.g0, self.sub_domains, 0)
        bc0 = [DirichletBC(V, g0, self.contour) for g0 in self.g0]

        # create function for inlet and outlet BC
        if self.waveform == "pulse":
            self.g1 = self.uExpr(t=t, **get_pulse_input_function(V, self.z_index, self.factor1, self.A, self.HR_inv, self.HR, self.b, self.f1))
            self.g2 = self.uExpr(t=t, **get_pulse_input_function(V, self.z_index, self.factor2, self.A, self.HR_inv, self.HR, self.b, self.f1))

        elif self.waveform =="sine":
            self.g1 = self.uExpr(t=t, **get_sine_input_function(V, self.z_index, self.HR, self.HR_inv, self.v_max1))
            self.g2 = self.uExpr(t=t, **get_sine_input_function(V, self.z_index, self.HR, self.HR_inv, self.v_max2))

        # Create inflow boundary condition for velocity on side 1 and 2
        bc1 = [DirichletBC(V, g1, self.top)    for g1 in self.g1]
        bc2 = [DirichletBC(V, g2, self.bottom) for g2 in self.g2]

        # Collect boundary conditions
        bcv = zip(bc1, bc0, bc2)
        bcp = [()]

        return bcv + bcp

    def initial_conditions(self, V, Q):
        u0 = self.uConstant((0.0, 0.0, 0.0))
        p0 = [Constant(0.0)]
        return u0 + p0

    def functional(self, t, u, p):

        print "time ", t
        v_max = max(_u.vector().norm("linf") for _u in as_list(u))
        print "max value of u ", v_max
        print "max value of p ", p.vector().norm("linf")

        f0 = self.flux(0,u)
        f1 = self.flux(1,u)
        f2 = self.flux(2,u)
        print "flux through A0 ", f0
        print "flux through A1 ", f1
        print "flux through A2 ", f2
        print "CFL = ", v_max * self.dt / self.h
        p_max = p.vector().norm("linf")
#       pressure_at_peak_v = p.vector()[self.pressure_index]

        return p_max

    def reference(self, t):
        # FIXME: Inset suitable reference value here
        return 0.0

    def flux(self, i, u):
        n = tetrahedron.n
        A = sum(u[d]*n[d] for d in range(3))*ds(i) # dot(u,n)*ds(i)
        return assemble(A, exterior_facet_domains=self.sub_domains, mesh=self.mesh)

    def area(self, i):
        f = Constant(1)
        A = f*ds(i)
        return assemble(A, exterior_facet_domains=self.sub_domains, mesh=self.mesh)

    def __str__(self):
        return "Chiari test problem"
