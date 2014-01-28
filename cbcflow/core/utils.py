

# --- Logging ---

from ..dol import MPI, warning, log, compile_extension_module
master = MPI.process_number() == 0

import os
from time import time

def cbcflow_warning(msg):
    if master:
        warning(msg)

def cbcflow_print(msg):
    if master:
        print msg

def cbcflow_log(level, msg):
    if master:
        log(level, msg)


def safe_mkdir(dir):
    """Create directory without exceptions in parallel."""
    # Create directory
    if not os.path.isdir(dir):
        try:
            os.makedirs(dir)
        except Exception as e:
            # Allow race condition when multiple processes
            # work in same directory, ignore exception.
            pass

    # Wait for all processes to finish, hopefully somebody
    # managed to create the directory...
    MPI.barrier()

    # Warn if this failed
    if not os.path.isdir(dir):
        #warning("FAILED TO CREATE DIRECTORY %s" % (dir,))
        Exception("FAILED TO CREATE DIRECTORY %s" % (dir,))
    

def timeit(t0=None, msg=None):
    if t0 is None:
        return time()
    else:
        t = time() - t0
        cbcflow_print("%s: %g" % (msg, t))
        return t

class Timer:
    def __init__(self, enabled):
        self._enabled = enabled
        self._timer = time()

    def completed(self, msg):
        t = time()
        ms = (t - self._timer)*1000
        if self._enabled:
            cbcflow_print("%10.0f ms: %s" % (ms, msg))
        self._timer = t


# --- System inspection ---

from os import getpid
from commands import getoutput
def get_memory_usage():
    mypid = getpid()
    mymemory = getoutput("ps -o rss %s" % mypid).split()[1]
    return mymemory


# --- Parallel hacks on top of dolfin ---

def parallel_eval(func, point, gather=True):
    """Parallel-safe function evaluation"""
    if gather:
        func.update()
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


# --- Network mesh retrieval ---

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

def retrieve(filename, urlbase='http://simula.no/~jobh/cbcflow'):
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


# --- String formatting ---

def time_to_string(t):
    t = int(t)
    s = ""
    hours = t/3600
    minutes = t/60-hours*60
    seconds = t-minutes*60-hours*3600

    if hours > 0: s += "%dh" %hours
    if minutes > 0: s += " %2dm" %minutes
    s += " %2ds" %seconds

    return s


# --- Type stuff ---

def as_list(u):
    "Return a list of objects."
    if isinstance(u, (list, tuple)):
        return u
    else:
        return [u]

def as_object(u):
    "Return a single object if possible, else a list."
    if not isinstance(u, (list, tuple)) or len(u) > 1:
        return u
    else:
        return u[0]


# --- Function space type stuff ---

def as_scalar_spaces(V):
    "Return a list of scalar (sub-)spaces consistent with V that can be used to apply DirichletBCs to components."
    d = V.cell().d
    if V.num_sub_spaces() == 0:
        return [V]*d
    else:
        return [V.sub(i) for i in xrange(d)]

def as_scalar_space(V):
    "Return a scalar (sub-)space consistent with V that can be used to construct a new Function."
    if V.num_sub_spaces() == 0:
        return V
    else:
        return V.sub(0).collapse()


# --- Maths ---
from ufl import grad, Identity
def epsilon(u):
    "Return symmetric gradient."
    return 0.5*(grad(u) + grad(u).T)

def sigma(u, p, mu):
    "Return stress tensor."
    return 2*mu*epsilon(u) - p*Identity(u.cell().d)


# --- Common solver parameters ---

maxiter = default_maxiter = 200
tolerance = default_tolerance = 1e-4


# --- Solution inspection ---

def has_converged(r, iter, method, maxiter=default_maxiter, tolerance=default_tolerance):
    "Check if solution has converged."
    cbcflow_print("Residual = %.3g" % r)
    if r < tolerance:
        cbcflow_print("%s iteration converged in %d iteration(s)." % (method, iter + 1))
        return True
    elif iter == maxiter - 1:
        raise RuntimeError("%s iteration did not converge." % method)
    return False

def is_periodic(bcs): # FIXME: Should we just remove this? Currently broken.
    "Check if boundary conditions are periodic."
    return False # FIXME: all(isinstance(bc, PeriodicBC) for bc in bcs)


# --- I/O stuff ---
def hdf5_link(hdf5filename, link_from, link_to):
    "Create internal link in hdf5 file"
    cpp_code = '''
    #include <hdf5.h>
    void link_dataset(const std::string hdf5_filename,
                              const std::string link_from,
                              const std::string link_to)
    {
        hid_t hdf5_file_id = HDF5Interface::open_file(hdf5_filename, "a", true);

        herr_t status;
        H5Lcreate_hard(hdf5_file_id, link_from.c_str(), H5L_SAME_LOC,
                       link_to.c_str(), H5P_DEFAULT, H5P_DEFAULT);
        dolfin_assert(status != HDF5_FAIL);
        HDF5Interface::close_file(hdf5_file_id);

    }
    '''

    cpp_module = compile_extension_module(cpp_code, additional_system_headers=["dolfin/io/HDF5Interface.h"])
    cpp_module.link_dataset(hdf5filename, link_from, link_to)