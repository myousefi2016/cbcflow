__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-04-03"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from dolfin import *

from time import time
import os
from os import getpid
from commands import getoutput

# Common solver parameters
maxiter = default_maxiter = 200
tolerance = default_tolerance = 1e-4

master = MPI.process_number() == 0

class SolverBase:
    "Base class for all solvers."

    def __init__(self, options):

        # Store options
        self.options = options

        # Reset some solver variables
        self._time = None
        self._cputime = 0.0
        self._timestep = 0

        # Reset files for storing solution
        self._ufiles = None
        self._pfile  = None

        # Reset storage for functional values and errors
        self._t = []
        self._M = []
        self._m = []
        self._e = []

        self._timer = None

    def getMyMemoryUsage(self):
        mypid = getpid()
        mymemory = getoutput("ps -o rss %s" % mypid).split()[1]
        return mymemory

    def timer(self, msg):
        if self.options['timer'] and master:
            print "%10.0f ms: %s"%((time()-self._timer)*1000, msg)
            self._timer = time()

    def start_timing(self):
        """Start timing, will be paused automatically during update
        and stopped when the end-time is reached."""
        self._time = time()
        self._timer = time()

    def solve(self, problem, dt, plot_solution=True):
        "Solve problem"
        raise NotImplementedError

    def prefix(self, problem):
        "Return file prefix for output files"
        p = problem.__module__.split(".")[-1].lower()
        s = self.__module__.split(".")[-1].lower()
        return problem.output_location + p + "_" + s

    def desegregate(self, u):
        if not isinstance(u, (list, tuple)):
            return u

        assert MPI.num_processes() == 1
        V = u[0].function_space()
        W = VectorFunctionSpace(V.mesh(), V.ufl_element().family(), V.ufl_element().degree())
        f = Function(W)
        fv = f.vector()
        for i in range(len(u)):
            uv = u[i].vector()
            # FIXME: Assumptions about ordering of vector-function entries
            fv[i*len(uv):(i+1)*len(uv)] = uv
        return f

    def update(self, problem, t, u, p):
        "Update problem at time t"

        casedir = os.path.join("results", self.options["casename"])

        if self.options['segregated']:
           s = max(ui.vector().norm('linf') for ui in u)
        else:
           s = u.vector().norm('linf')
        if s > 5 * getattr(problem, 'U', float('inf')):
            warning("A component in u is %.4g times characteristic velocity U"%round(s))
        if s > 1e10:
            raise RuntimeError("Runaway solution")

        # Add to accumulated CPU time
        timestep_cputime = time() - self._time
        self._cputime += timestep_cputime

        # Compute divergence
        if self.options["compute_divergence"]:
            check_divergence(u, p.function_space())

        # Update problem FIXME: Should this be called before problem.functional??
        problem.update_problem(t, self._list_or_function(u), p)

        # Evaluate functional and error
        m = problem.reference(t)
        M = problem.functional(t, self._list_or_function(u), p)
        if m is None:
            e = None
            if master:
                print "M = %g (missing reference value)" % M
        else:
            e = abs(M - m)
            if master:
                print "M = %g (reference %g), error = %g (maximum %g)" % (M, m, e, max([e] + self._e))

        # Store values
        self._t.append(t)
        self._M.append(M)
        self._m.append(m)
        self._e.append(e)

        # Save solution
        if self.options["save_solution"]:

            # Save velocity and pressure
            frequency = self.options["save_frequency"]

            if (self._timestep - 1) % frequency == 0:
                # Create files
                if self._ufiles is None:
                    if self.options['segregated']:
                        # added '_' to separate i from vtu numbering
                        self._ufiles = [File(os.path.join(casedir, "u%d_.pvd" % i))
                                        for i in range(len(u))]
                    else:
                        self._ufiles = File(os.path.join(casedir, "u.pvd"))
                if self._pfile is None:
                    self._pfile = File(os.path.join(casedir, "p.pvd"))

                # Write to files
                if self.options['segregated']:
                    for i, ui in enumerate(u):
                        self._ufiles[i] << ui
                else:
                    self._ufiles << u[0]
                self._pfile << p

        # Save solution at t = T
        if self.options["save_solution_at_t=T"]:
            if t >= problem.T:
                # Create files
                if self._ufiles is None:
                    self._ufiles = [File(os.path.join(casedir, "u%d_at_end.pvd" % i))
                                    for i in range(len(u))]
                else: 
                    self._ufiles = File(os.path.join(casedir, "u_at_end.pvd"))
                if self._pfile is None:
                    self._pfile = File(os.path.join(casedir, "p_at_end.pvd"))

                # Write to files
                if self.options['segregated']:
                    for i, ui in enumerate(u):
                        self._ufiles[i] << ui
                else:
                    self._ufiles << u[0]
                self._pfile << p

        # Save vectors in xml format
        if self.options["save_xml"]:
            timestr = "at_t%d_%.5e" % (self._timestep, t)
            if self.options['segregated']:
                for i, ui in enumerate(u):
                    file = File(os.path.join(casedir, "u%d_%s.xml.gz" % (i, timestr)))
                    file << ui.vector()
            else:
                file = File(os.path.join(casedir, "u_%s.xml.gz" % (timestr,)))
                file << u[0].vector()
            file = File(os.path.join(casedir, "p_%s.xml.gz" % (timestr,)))
            file << p.vector()

        # Plot solution
        if self.options["plot_solution"]:
            # Plot velocity and pressure
            plot(self.desegregate(u), title="Velocity", rescale=True)
            plot(p, title="Pressure", rescale=True)

        # Check memory usage
        if self.options["check_mem_usage"]:
            if (self._timestep - 1) % self.options["check_frequency"] == 0:
                print 'Memory usage is:' , self.getMyMemoryUsage()

        # Print progress
        if master:
            ss = self._cputime * (problem.T/t-1)
            hh, ss = divmod(ss, 60*60)
            mm, ss = divmod(ss, 60)
            print
            s = "Time step %d finished in %.2f seconds, %.1f%% done (t=%.3g, T=%g; %02d:%02d:%02d remaining)." \
                % (self._timestep, timestep_cputime, 100.0*(t / problem.T), t, problem.T, hh, mm, ss)
            print s
            print "-"*len(s)

        # Increase time step and record current time
        self._timestep += 1
        self._time = time()

    def eval(self):
        "Return last functional value and maximum error in functional value on [0, T]"

        # Plot values
        if self.options["plot_functional"]:
            from pylab import plot, xlabel, ylabel, grid, show
            plot(self._t, self._M)
            xlabel('t')
            ylabel('Functional')
            grid(True)
            show()

        # Return value
        if self._e[0] is None:
            return self._M[-1], None
        else:
            return self._M[-1], max([0.0] + self._e)

    def cputime(self):
        "Return accumulated CPU time."
        return self._cputime

    def _list_or_function(self, u):
        "Return a single function if possible, else a list."
        if self.options['segregated'] or not isinstance(u, (list, tuple)):
            return u
        else:
            return u[0]

def epsilon(u):
    "Return symmetric gradient."
    return 0.5*(grad(u) + grad(u).T)

def sigma(u, p, nu):
    "Return stress tensor."
    return 2*nu*epsilon(u) - p*Identity(u.cell().d)

def is_periodic(bcs):
    "Check if boundary conditions are periodic."
    return all(isinstance(bc, PeriodicBC) for bc in bcs)

def has_converged(r, iter, method, maxiter=default_maxiter, tolerance=default_tolerance):
    "Check if solution has converged."
    if master: print "Residual = ", r
    if r < tolerance:
        if master: print "%s iteration converged in %d iteration(s)." % (method, iter + 1)
        return True
    elif iter == maxiter - 1:
        raise RuntimeError, "%s iteration did not converge." % method
    return False

def check_divergence(u, Q):
    "Check divergence of velocity."

    # Compute L2 norm of divergence
    if master: print "||div u||_L2 =", norm(u, "Hdiv0")

    # Compute projection of div u into Q_0
    pdivu = project(div(u), Q)
    zero = Constant(Q.mesh(), 0.0)
    bc = DirichletBC(Q, zero, DomainBoundary())
    bc.apply(pdivu.vector())

    # Compute "weak" L2 norm of divergence
    if master: print "||div u||_w  =", sqrt(abs(assemble(pdivu*div(u)*dx, mesh=Q.mesh())))
