from __future__ import division
__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-03-19"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kent-Andre Mardal, 2008.

from dolfin import *
from math import *
import os

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

    def pressure_bc(self, Q):
        if master:
            warning("Using default pressure 0, please set pressure bc in boundary_conditions()")
        return Constant(0)

    def uConstant(self, values):
        if isinstance(values, tuple) and self.options['segregated']:
            return [Constant(v) for v in values]
        else:
            return [Constant(values)]

    def uExpr(self, cppcode, **kwargs):
        if self.options['segregated'] and isinstance(cppcode, (tuple, list)):
            return [Expression(e, **kwargs) for e in cppcode]
        else:
            return [Expression(cppcode, **kwargs)]

    def eval(self, func, point, gather=True):
        """Parallel-safe function evaluation"""
        if gather:
            if hasattr(func, 'update'):
                func.update() # dolfin dev
            else:
                func.gather() # dolfin 1.0
        if len(func.shape())==1:
            M = [0]*func.shape()[0]
        else:
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

    def retrieve(self, filename, urlbase='http://simula.no/~jobh/headflow'):
        if not filename.endswith(".gz"):
            # Enforcing .gz extension is a quick fix to avoid trouble when
            # httpserver serves .gz file without extension, which is then
            # unreadable for dolfin.
            filename += ".gz"
        if master and not os.path.exists(filename):
            url = urlbase+'/'+filename
            warning('%s not found, fetching from %s'%(filename,url))

            targetdir = os.path.abspath(filename[:filename.rfind('/')])
            log_level = get_log_level()
            set_log_level(PROGRESS)
            progress = [Progress(filename.split('/')[-1])]
            def reporter(numblocks, blocksize, totalsize):
                progress[0] += numblocks*blocksize / totalsize

            if not os.path.isdir(targetdir):
                os.makedirs(targetdir)
            try:
                DataURLOpener(url, filename).retrieve(reporter)
            except:
                if os.path.exists(filename):
                    os.remove(filename)
                raise

            del progress[0]
            set_log_level(log_level)

        MPI.barrier()
        return filename

def as_list(u):
    "Return a list of objects."
    if isinstance(u, (list, tuple)):
        return u
    else:
        return [u]

import urllib
class DataURLOpener(urllib.FancyURLopener):
    def __init__(self, url, filename):
        urllib.FancyURLopener.__init__(self)
        self.url = url
        self.filename = filename
    def retrieve(self, reporter=None, data=None):
        urllib.FancyURLopener.retrieve(self, self.url, self.filename, reporter, data)
    def http_error_default(self, url, fp, errcode, errmsg, headers):
        raise IOError(str(errcode)+" "+errmsg+", "+self.url)
