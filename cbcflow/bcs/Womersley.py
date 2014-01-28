
from ..dol import *

import numpy as np

from scipy.interpolate import UnivariateSpline
from scipy.integrate import simps
from scipy.special import jn

from itertools import izip
from .bcutils import compute_boundary_geometry_acrn, compute_transient_scale_value, x_to_r2

def fourier_coefficients(x, y, T, N=25):
    '''From x-array and y-spline and period T, calculate N complex Fourier coefficients.'''
    omega = 2*pi/T
    ck = []
    ck.append(1/T*simps(y(x), x))
    for n in range(1,N):
        c = 1/T*simps(y(x)*np.exp(-1j*n*omega*x), x)

        # Clamp almost zero real and imag components to zero
        if 1:
            cr = c.real
            ci = c.imag
            if abs(cr) < 1e-14: cr = 0.0
            if abs(ci) < 1e-14: ci = 0.0
            c = cr + ci*1j

        ck.append(2*c)
    return ck

class WomersleyComponent2(Expression):
    # Subclassing the expression class restricts the number of arguments, args is therefore a dict of arguments.
    def __init__(self, args): # TODO: Document args properly
        Expression.__init__(self)

        # Spatial args
        self.radius = args["radius"]
        self.center = args["center"]
        self.normal = args["normal"]
        self.normal_component = args["normal_component"]

        # Temporal args
        self.period = args["period"]
        self.Qn = args["Qn"]

        # Physical args
        self.nu = args["nu"]

        # Internal state
        self._t = None
        self._scale_value = 1.0

        # Precomputation
        self._values = np.zeros(10000)
        self._ys = np.linspace(0.0, 1.0, len(self._values))
        self._precompute_bessel_functions()
        self._precompute_y_coeffs()

    def _precompute_bessel_functions(self):
        '''Calculate the Bessel functions of the Womersley profile'''
        self.omega = 2 * pi / self.period
        self.N = len(self.Qn)
        self.ns = np.arange(1, self.N)

        # Allocate for 0...N-1
        alpha = np.zeros(self.N, dtype=np.complex)
        self.beta = np.zeros(self.N, dtype=np.complex)
        self.jn0_betas = np.zeros(self.N, dtype=np.complex)
        self.jn1_betas = np.zeros(self.N, dtype=np.complex)

        # Compute vectorized for 1...N-1 (keeping element 0 in arrays to make indexing work out later)
        alpha[1:] = self.radius * np.sqrt(self.ns * (self.omega / self.nu))
        self.beta[1:] = alpha[1:] * np.sqrt(1j**3)
        self.jn0_betas[1:] = jn(0, self.beta[1:])
        self.jn1_betas[1:] = jn(1, self.beta[1:])

    def _precompute_y_coeffs(self):
        "Compute intermediate terms for womersley function."
        n = len(self._values)
        self._y_coeffs = np.zeros((n,self.N), dtype=np.complex)
        pir2 = pi * self.radius**2
        for i, y in enumerate(self._ys):
            self._y_coeffs[i,0] = (2*self.Qn[0]/pir2) * (1 - y**2)
            for n in self.ns:
                tmp1 = self.Qn[n] / pir2
                tmp2 = 1.0 - jn(0, self.beta[n]*y) / self.jn0_betas[n]
                tmp3 = 1.0 - 2.0*self.jn1_betas[n] / (self.beta[n]*self.jn0_betas[n])
                self._y_coeffs[i,n] = tmp1 * (tmp2 / tmp3)

    def set_t(self, t):
        # Compute time dependent coeffs once
        self._t = float(t) % self.period
        self._expnt = np.exp((self.omega * self._t * 1j) * self.ns)

        # Compute values for this time for each y
        n = len(self._values)
        for i in xrange(n):
            # Multiply complex coefficients for x with complex exponential functions in time
            wom = (self._y_coeffs[i,0] + np.dot(self._y_coeffs[i,1:], self._expnt)).real

            # Scale by negative normal direction and scale_value
            self._values[i] = -self.normal_component * self._scale_value * wom

    def eval(self, value, x):
        y = np.sqrt(x_to_r2(x, self.center, self.normal)) / self.radius
        nm = len(self._values) - 1
        yi = max(0, min(nm, int(round(y*nm))))
        #print "w(%.2e, %.2e, %.2e) = %.2e // y = %.2e" % (x[0], x[1], x[2], val, y)
        value[0] = self._values[yi]

class WomersleyComponent1(Expression):
    # Subclassing the expression class restricts the number of arguments, args is therefore a dict of arguments.
    def __init__(self, args): # TODO: Document args properly
        Expression.__init__(self)

        # Spatial args
        self.radius = args["radius"]
        self.center = args["center"]
        self.normal = args["normal"]
        self.normal_component = args["normal_component"]

        # Temporal args
        self.period = args["period"]
        if "Q" in args:
            assert "V" not in args, "Cannot provide both Q and V!"
            self.Qn = args["Q"]
            self.N = len(self.Qn)
        elif "V" in args:
            self.Vn = args["V"]
            self.N = len(self.Vn)
        else:
            error("Invalid transient data type, missing argument 'Q' or 'V'.")

        # Physical args
        self.nu = args["nu"]

        # Internal state
        self.t = None
        self.scale_value = 1.0

        # Precomputation
        self._precompute_bessel_functions()
        self._all_r_dependent_coeffs = {}

    def _precompute_bessel_functions(self):
        '''Calculate the Bessel functions of the Womersley profile'''
        self.omega = 2 * pi / self.period
        self.ns = np.arange(1, self.N)

        # Allocate for 0...N-1
        alpha = np.zeros(self.N, dtype=np.complex)
        self.beta = np.zeros(self.N, dtype=np.complex)
        self.jn0_betas = np.zeros(self.N, dtype=np.complex)
        self.jn1_betas = np.zeros(self.N, dtype=np.complex)

        # Compute vectorized for 1...N-1 (keeping element 0 in arrays to make indexing work out later)
        alpha[1:] = self.radius * np.sqrt(self.ns * (self.omega / self.nu))
        self.beta[1:] = alpha[1:] * np.sqrt(1j**3)
        self.jn0_betas[1:] = jn(0, self.beta[1:])
        self.jn1_betas[1:] = jn(1, self.beta[1:])

    def _precompute_r_dependent_coeffs(self, y):
        pir2 = pi * self.radius**2
        # Compute intermediate terms for womersley function
        r_dependent_coeffs = np.zeros(self.N, dtype=np.complex)
        if hasattr(self, 'Vn'):
            #r_dependent_coeffs[0] = (self.Vn[0]/2.0) * (1 - y**2)
            r_dependent_coeffs[0] = self.Vn[0] * (1 - y**2)
            for n in self.ns:
                r_dependent_coeffs[n] = self.Vn[n] * (self.jn0_betas[n] - jn(0, self.beta[n]*y)) / (self.jn0_betas[n] - 1.0)
        elif hasattr(self, 'Qn'):
            r_dependent_coeffs[0] = (2*self.Qn[0]/pir2) * (1 - y**2)
            for n in self.ns:
                bn = self.beta[n]
                j0bn = self.jn0_betas[n]
                j1bn = self.jn1_betas[n]
                r_dependent_coeffs[n] = (self.Qn[n] / pir2) * (j0bn - jn(0, bn*y)) / (j0bn - (2.0/bn)*j1bn)
        else:
            error("Missing Vn or Qn!")
        return r_dependent_coeffs

    def _get_r_dependent_coeffs(self, y):
        "Look for cached womersley coeffs."
        key = y
        r_dependent_coeffs = self._all_r_dependent_coeffs.get(key)
        if r_dependent_coeffs is None:
            # Cache miss! Compute coeffs for this coordinate the first time.
            r_dependent_coeffs = self._precompute_r_dependent_coeffs(y)
            self._all_r_dependent_coeffs[key] = r_dependent_coeffs
        return r_dependent_coeffs

    def set_t(self, t):
        self.t = float(t) % self.period
        self._expnt = np.exp((self.omega * self.t * 1j) * self.ns)

    def eval(self, value, x):
        # Compute or get cached complex coefficients that only depend on r
        y = np.sqrt(x_to_r2(x, self.center, self.normal)) / self.radius
        coeffs = self._get_r_dependent_coeffs(y)

        # Multiply complex coefficients for x with complex exponential functions in time
        wom = (coeffs[0] + np.dot(coeffs[1:], self._expnt)).real

        # Scale by negative normal direction and scale_value
        value[0] = -self.normal_component * self.scale_value * wom
        #print "w(%.2e, %.2e, %.2e) = %.2e" % (x[0], x[1], x[2], value[0])


def make_womersley_bcs(coeffs, mesh, indicator, nu, scale_to=None, facet_domains=None,
                       coeffstype="Q", num_fourier_coefficients=25):
    """Generate a list of expressions for the components of a Womersley profile."""
    assert(isinstance(mesh, Mesh))

    # TODO: Always require facet_domains
    if facet_domains is None:
        dim = mesh.geometry().dim()
        facet_domains = MeshFunction("size_t", mesh, dim-1, mesh.domains())

    # Compute boundary geometry
    area, center, radius, normal = compute_boundary_geometry_acrn(mesh, indicator, facet_domains)

    # Compute transient profile as interpolation of given coefficients
    x,y = zip(*coeffs)
    x = np.array(x)
    y = np.array(y)
    period = max(x)
    transient_profile = UnivariateSpline(x, y, s=0, k=1)

    # Compute fourier coefficients of transient profile
    timedisc = np.linspace(0, period, 1001)
    Cn = fourier_coefficients(timedisc, transient_profile, period, num_fourier_coefficients)
    if 0: # FIXME: Move this code into a unit test of fourier_coefficients:
        print "*"*80
        print "Cn =", Cn
        print "Reconstructing transient profile from Cn:"
        for t in np.linspace(0.0, period, 100):
            evaluated = transient_profile(t)
            reconstructed = np.dot(Cn, np.exp(1j*(2*pi/period)*np.arange(len(Cn))*t)).real
            print "%.2e  %.2e  %.2e  %.2e" % (evaluated, reconstructed, (evaluated - reconstructed),
                                              (evaluated - reconstructed) * 2 / abs(evaluated + reconstructed))
        print "*"*80

    if 0:
        print "*"*80
        print "In womersley:"
        print 'r', radius
        print 'c', center
        print 'n', normal
        print 'om', period
        print 'nu=', nu
        print "*"*80

    # Create Expressions for each direction
    expressions = []
    for ncomp in normal:
        args = {
            "radius": radius,
            "center": center,
            "normal": normal,
            "normal_component": ncomp,
            "period": period,
            "nu": nu,
            }
        args[coeffstype] = Cn
        expressions.append(WomersleyComponent1(args))

    # Apply scaling w.r.t. peak transient profile (FIXME: This is unclear!)
    if scale_to is not None:
        scale_factor = compute_transient_scale_value(expressions, period,
                                                     mesh, facet_domains, indicator,
                                                     scale_to)
        for e in expressions:
            e.scale_value = scale_factor

    return expressions

class Womersley(list):
    def __init__(self, coeffs, mesh, indicator, nu, scale_to=None, facet_domains=None):
        print "Deprecation warning: use make_womersley_bcs instead of Womersley class." # FIXME: Remove class
        self.extend(make_womersley_bcs(coeffs, mesh, indicator, nu, scale_to, facet_domains))