
# Always import the entire dolfin namespace
from dolfin import *

# If dolfin_adjoint has already been imported, overwrite the
# dolfin namespace with the overloaded dolfin_adjoint symbols
import sys

has_dolfin_adjoint = ('dolfin_adjoint' in sys.modules)
if has_dolfin_adjoint:
    from dolfin_adjoint import *
else:
    def adj_start_timestep(*args):
        pass
    def adj_inc_timestep(*args, **kwargs):
        pass

    # TODO: Do we need something like this? Keeping this code for a while just in case:
    #def project(*args, **kwargs):
    #    if "name" in kwargs: del kwargs["name"]
    #    if "annotate" in kwargs: del kwargs["annotate"]
    #    return dolfin.project(*args, **kwargs)
    #def assemble(*args, **kwargs):
    #    if "name" in kwargs: del kwargs["name"]
    #    if "annotate" in kwargs: del kwargs["annotate"]
    #    return dolfin.assemble(*args, **kwargs)

# This is a trick to handle automatic timestep annotation
def Time(t0=0.0):
    t = Constant(t0, name="TIME")
    t._prev_value = t0
    t._assigned_to = False
    return t

def assign_time(t, tvalue):
    if t._assigned_to:
        # Annotate previous timestep is done
        adj_inc_timestep(t._prev_value)
    else:
        # Annotate the beginning of time
        t._assigned_to = True
        adj_start_timestep(t._prev_value)
    # Update time constant to reflect modern times
    t.assign(tvalue)
    t._prev_value = tvalue

def finalize_time(t):
    # Make sure we have annotated the beginning of time
    if not t._assigned_to:
        t._assigned_to = True
        adj_start_timestep(t._prev_value)
    # Annotate the end-time is here
    adj_inc_timestep(t._prev_value, finished=True)
    # Time constant needs no updating anymore
    t._prev_value = None