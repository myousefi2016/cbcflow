
from ..dol import *

import numpy as np
from scipy.interpolate import UnivariateSpline
from itertools import izip
from .bcutils import compute_boundary_geometry_acrn, compute_transient_scale_value, x_to_r2

class PoiseuilleComponent(Expression):
    # Subclassing the expression class restricts the number of arguments, args is therefore a dict of arguments.
    def __init__(self, args): # TODO: Document args properly
        Expression.__init__(self)

        # Spatial args
        self.radius = args["radius"]
        self.center = args["center"]
        self.normal = args["normal"]
        self.normal_component = args["normal_component"]

        if 0:
            print "RADIUS:", self.radius
            print "CENTER:", self.center
            print "NORMAL:", self.normal_component

        # Temporal args
        self.period = args["period"]
        # FIXME: Remove this temporary hack, here to allow external problems to be updated smoothly:
        if "velocity_profile" in args:
            print "NB! 'velocity_profile' is deprecated, renamed to transient_profile!"
            self.transient_profile = args["velocity_profile"]
        else:
            self.transient_profile = args["transient_profile"]

        # Internal state
        self.t = 0.0
        self.scale_value = 1.0

    def set_t(self, t):
        self.t = float(t) % self.period
        self._tp = 2.0 * self.transient_profile(self.t) / (pi * self.radius**2)

    def eval(self, value, x):
        # Compute radial coordinates
        #r2 = sum((xi-ci)**2 for xi,ci in izip(x,self.center))
        r2 = x_to_r2(x, self.center, self.normal)
        y2 = r2 / self.radius**2

        # Compute scalar velocity profile value in flow direction
        velocity_profile = (1 - y2)

        # Scale by negative normal direction, scale_value, and transient profile
        val = -self.normal_component * self.scale_value * velocity_profile * self._tp

        # Output final value
        value[0] = val


def make_poiseuille_bcs(coeffs, mesh, indicator, scale_to=None, facet_domains=None):
    """Generate a list of expressions for the components of a Poiseuille profile."""
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

    if 0:
        print "*"*80
        print "In poiseuille:"
        print 'r', radius
        print 'c', center
        print 'n', normal
        print 'om', period
        print 'Q(0.4)', transient_profile(0.4)
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
            "transient_profile": transient_profile,
            }
        expressions.append(PoiseuilleComponent(args))

    # Apply scaling w.r.t. peak transient profile (FIXME: This is unclear!)
    if scale_to is not None:
        scale_factor = compute_transient_scale_value(expressions, period,
                                                     mesh, facet_domains, indicator,
                                                     scale_to)
        for e in expressions:
            e.scale_value = scale_factor

    return expressions

class Poiseuille(list):
    def __init__(self, coeffs, mesh, indicator, scale_to=None, facet_domains=None):
        print "Deprecation warning: use make_poiseuille_bcs instead of Pouseille class." # FIXME: Remove class
        self.extend(make_poiseuille_bcs(coeffs, mesh, indicator, scale_to, facet_domains))
