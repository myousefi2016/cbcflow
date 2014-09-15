# Copyright (C) 2010-2014 Simula Research Laboratory
#
# This file is part of CBCFLOW.
#
# CBCFLOW is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CBCFLOW is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CBCFLOW. If not, see <http://www.gnu.org/licenses/>.
'''
import os
from time import time

from cbcflow.dol import MPI, warning, log, compile_extension_module

def on_master_process():
    return MPI.process_number() == 0

def in_serial():
    return MPI.num_processes() == 1


# --- Logging ---

def cbc_warning(msg):
    if on_master_process():
        warning(msg)

def cbc_print(msg):
    if on_master_process():
        print msg

def cbc_log(level, msg):
    if on_master_process():
        log(level, msg)


def safe_mkdir(dir):
    """Create directory without exceptions in parallel."""
    # Create directory
    if not os.path.isdir(dir):
        try:
            os.makedirs(dir)
        except:
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
        cbc_print("%s: %g" % (msg, t))
        return t

class Timer:
    def __init__(self, frequency=0):
        self._frequency = frequency
        self._timer = time()
        self._timings = {}
        self._keys = []
        self._N = 0

    def completed(self, key, summables={}):
        if self._frequency == 0:
            return
        
        if key not in self._timings:
            self._keys.append(key)
            self._timings[key] = [0,0, {}]
        
        t = time()
        ms = (t - self._timer)*1000
        self._timings[key][0] += ms
        self._timings[key][1] += 1
        
        for k,v in summables.items():
            if k not in self._timings[key][2]:
                self._timings[key][2][k] = 0
            self._timings[key][2][k] += v
        
        if self._frequency == 1:
            s = "%10.0f ms: %s" % (ms, key)
            ss = []
            #if summables != {}:
            #    s += "  ("
            for k, v in summables.items():
                #s += "%s: %s, " %(k,v)
                ss.append("%s=%s" %(k,v))
            if len(ss) > 0:
                ss = "  ("+", ".join(ss)+")"
            else:
                ss = ""
            s += ss

            cbc_print(s)

        self._timer = time()
    
    def _print_summary(self):
        cbc_print("Timings summary: ")
        
        for key in self._keys:
            tot = self._timings[key][0]
            N = self._timings[key][1]
            avg = int(1.0*tot/N)
            
            s = "%10.0f ms (avg: %8.0f ms, N: %5d): %s" %(tot, avg, N, key)
            
            summables = self._timings[key][2]
            ss = []
            #if summables != {}:
            #    s += "("
            for k, tot in summables.items():
                avg = int(1.0*tot/N)
                #s += "%s: %s (avg: %s), " %(k,tot,avg)
                ss.append("%s=%s (avg: %s)" %(k,tot,avg))
            #if summables != {}:
            #    s += ")"
            if len(ss) > 0:
                ss = "  ("+", ".join(ss)+")"
            else:
                ss = ""
            s += ss
            cbc_print(s)
    
    def _reset(self):
        self._timings = {}
        self._keys = []
        self._N = 0
        
    def increment(self):
        self._N += 1
        if self._frequency > 1 and self._N % self._frequency == 0:
            self._print_summary()
            self._reset()
    

# --- System inspection ---

from os import getpid
from commands import getoutput
def get_memory_usage():
    """Return memory usage in MB"""
    try:
        from fenicstools import getMemoryUsage
        return getMemoryUsage()
    except:
        cbc_warning("Unable to load fenicstools to check memory usage. Falling back to unsafe memory check.")
        mypid = getpid()
        mymemory = getoutput("ps -o rss %s" % mypid).split()[1]
        return int(mymemory)/1024


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
    if on_master_process() and N > 1:
        warning("%d processors returned function value, which is unexpected (but probably ok)"%N)
    if hasattr(M, '__iter__'):
        for i in range(len(M)):
            M[i] = MPI.sum(M[i])/N
    else:
        M = MPI.sum(M)/N
    return M
'''
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


"""
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
"""

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
    cbc_print("Residual = %.3g" % r)
    if r < tolerance:
        cbc_print("%s iteration converged in %d iteration(s)." % (method, iter + 1))
        return True
    elif iter == maxiter - 1:
        raise RuntimeError("%s iteration did not converge." % method)
    return False

def is_periodic(bcs): # FIXME: Should we just remove this? Currently broken.
    "Check if boundary conditions are periodic."
    return False # FIXME: all(isinstance(bc, PeriodicBC) for bc in bcs)
