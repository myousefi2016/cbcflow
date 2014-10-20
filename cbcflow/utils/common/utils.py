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
    return 2*mu*epsilon(u) - p*Identity(u.cell().geometric_dimension())


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

from dolfin import Function, FunctionAssigner
from cbcflow.utils.core.spaces import NSSpacePoolMixed, NSSpacePoolSegregated
class PressureConverter():
    def __call__(self, p, spaces):
        if not isinstance(p, Function):
            if not hasattr(self, "_p"):
                self._p = Function(spaces.Q)
                assert isinstance(spaces, NSSpacePoolMixed)
                self._assigner = FunctionAssigner(spaces.Q, spaces.W.sub(1))

            # Hack: p is a Indexed(Coefficient()),
            # get the underlying mixed function
            w = p.operands()[0]
            self._assigner.assign(self._p, w.sub(1))

            p = self._p

        assert isinstance(p, Function)
        return p
    
class VelocityConverter():
    def __call__(self, u, spaces):
        if not isinstance(u, Function):
            d = spaces.d
            if not hasattr(self, "_u"):
                self._u = Function(spaces.V)

                if isinstance(spaces, NSSpacePoolMixed):
                    self._assigner = FunctionAssigner(spaces.V, spaces.W.sub(0))
                elif isinstance(spaces, NSSpacePoolSegregated):
                    self._assigner = FunctionAssigner(spaces.V, [spaces.U]*d)
                else:
                    error("It doesnt make sense to create a function assigner for a split space.")

            if isinstance(spaces, NSSpacePoolMixed):
                # Hack: u is a ListTensor([Indexed(Coefficient()),...]),
                # get the underlying mixed function
                w = u.operands()[0].operands()[0]
                assert w.shape() == (d+1,)
                us = w.sub(0)

            elif isinstance(spaces, NSSpacePoolSegregated):
                us = [u[i] for i in range(d)]

            self._assigner.assign(self._u, us)
            u = self._u

        assert isinstance(u, Function)
        return u

