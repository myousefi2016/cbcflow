__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-03-19"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kent-Andre Mardal, 2008.

from dolfin import *
from numpy import linspace
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

    # FIXME: Move to solverbase.py and rename select_timestep

    def timestep(self, problem):
        "Return time step and number of time steps for problem."

        # FIXME: This looks very complex, should be cleaned up

        T  = problem.T
        U  = problem.U
        nu = problem.nu
        h  = MPI.min(self.mesh.hmin())
        "Return time step and number of time steps for problem. Used for debugging / compilation only"
        if self.options["max_steps"] is not None:
            dt =  0.25*h**2 / (U*(nu + h*U))
            n  = self.options["max_steps"]
            T  = n*dt
            t_range = linspace(0, T, n + 1)[1:]
        else:
            # FIXME: This sequence of ifs make no sense. Clean up...

            if self.options["dt"]:
                dt = self.options["dt"]
                n = int(T / dt + 0.5)
                dt = T / n
                if master: print "Using used supplied dt"

            elif self.options["dt_division"] != 0 and getattr(problem, "dt", 0) > 0:
                dt = problem.dt / int(sqrt(2)**self.options["dt_division"])
                n  = int(T / dt + 1.0)
                dt = T / n
                if master: print 'Using problem.dt and time step refinements'

            # Use time step specified in problem if available
            elif getattr(problem, "dt", 0) > 0:
                dt = problem.dt
                n  = int(T / dt)
                if master: print 'Using problem.dt'

            # Otherwise, base time step on mesh size
            elif self.options["dt_division"] != 0:
                dt = 0.25*h**2 / (U*(nu + h*U))
                dt /= int(sqrt(2)**self.options["dt_division"])
                n  = int(T / dt + 1.0)
                dt = T / n
                if master: print 'Computing time step according to stability criteria and time step refinements'

            # Otherwise, base time step on mesh size
            else:
#                dt =  0.25*h**2 / (U*(nu + h*U))
                dt =  0.2*(h / U)
                n  = int(T / dt + 1.0)
                dt = T / n
                if master: print 'Computing time step according to stability criteria'

        # Compute range
        t_range = linspace(0,T,n+1)[1:] # FIXME: Comment out [1:] to run g2ref g2ref

        # Report time step
        if master:
            print " "
            print 'Number of timesteps:' , len(t_range)
            print 'Size of timestep:' , dt
            print " "

        return dt, t_range[0], t_range

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

    def uEval(self, func, component, point):
        if self.options['segregated']:
            return func[component](point)
        else:
            return func(point)[component]
