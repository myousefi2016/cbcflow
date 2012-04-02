__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-03-19"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kent-Andre Mardal, 2008.

from dolfin import *
from math import *

master = MPI.process_number() == 0

class ProblemBase:
    "Base class for all problems."

    def __init__(self, options):

        # Store options
        self.options = options

        self.output_location = ''

        # FIXME: Temporary while testing
        self.bcpsi = []

    def update_problem(self, t, u, p):
        "Update problem at time t"

        # Update state
        self.t = t
        self.u = u
        self.p = p

        # Call problem-specific update
        self.update(t, u, p)

    def update(self, t, u, p):
        "Problem-speficic update at time t"
        pass

    def functional(self, t, u, p):
        "Return value of functional of interest"
        return 0.0

    def reference(self, t):
        "Return reference value for functional"
        return None

    def tolerance(self, problem):
        "Return tolerance (used as local convergence criterion)."
        if str(problem) == 'Channel':
            return 1e-11
        elif str(problem) == 'Cylinder':
            return 1e-7
        else:
            return 1e-6

    def resistance(self, mesh, boundary_markers, mark, C, p0):
        if not hasattr(self, "u"):
            if master:
                print "self.u not initialized, assuming zero flux (resistance is %.3g)"%p0
            return Constant(p0)
        n = FacetNormal(mesh)
        flux = inner(self.u, n)*ds(mark)
        Q = assemble(flux, mesh=mesh, exterior_facet_domains=boundary_markers)
        R = C*Q + p0
        if master:
            print "Computed resistance over marker %d is %.3g, the flux is %.3g"%(mark, R, Q)
        return Constant(R)

    def pressure_bc(self, Q):
        if master:
            warning("Using default pressure 0, please set pressure bc in boundary_conditions()")
        return Constant(0)

    def uConstant(self, values):
        if isinstance(values, tuple) and self.options['segregated']:
            return [Constant(v) for v in values]
        else:
            return [Constant(values)]

    def uExpr(self, expr_strings, **kwargs):
        if self.options['segregated'] and isinstance(expr_strings, (tuple, list)):
            return [Expression(e, **kwargs) for e in expr_strings]
        else:
            return [Expression(expr_strings, **kwargs)]

    def eval(self, func, point, gather=True):
        """Parallel-safe function evaluation"""
        if gather:
            if hasattr(func, 'gather'):
                func.gather() # dolfin 1.0
            else:
                func.update() # dolfin dev
        M = 0
        try:
            M = func(point)
            N = MPI.sum(1) # Succeeding processors participate in the MPI collective here
        except RuntimeError:
            N = MPI.sum(0) # Failing processors participate in the MPI collective here
            if N == 0:
                raise      # All processors failed
        if master and N > 1:
            warning("%d processors returned function value, which is unexpected (but probably ok)"%N)
        if hasattr(M, '__iter__'):
            for i in range(len(M)):
                M[i] = MPI.sum(M[i])/N
        else:
            M = MPI.sum(M)/N
        return M

    def uEval(self, func, component, point):
        if self.options['segregated']:
            return self.eval(func[component], point)
        else:
            return self.eval(func, point)[component]
