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

from cbcflow.core.parameterized import Parameterized


class DependencyException(Exception):
    def __init__(self, fieldname=None, dependency=None, timestep=None, original_exception_msg=None):
        message = []
        if fieldname:
            message += ["Dependency/dependencies not found for field %s." % fieldname]
        if dependency:
            message += ["Dependency %s not functioning." % dependency]
        if timestep:
            message += ["Relative timestep is %d. Are you trying to calculate time-derivatives at t=0?" % (timestep)]
        if original_exception_msg:
            message += ["\nOriginal exception was: " + original_exception_msg]
        message = ' '.join(message)
        Exception.__init__(self, message)


# Fields available through pp.get(name) even though they have no Field class
builtin_fields = ("t", "timestep")

class PostProcessor(Parameterized):
    def __init__(self, params=None):
        Parameterized.__init__(self, params)

        # Storage of actual fields
        self._fields = {}

        # Representations of field dependencies
        self._sorted_fields_keys = [] # Topological ordering of field names
        self._dependencies = {} # Direct dependencies dep[name] = ((depname0, ts0), (depname1, ts1), ...)
        self._full_dependencies = {} # Indirect dependencies included
        for depname in builtin_fields:
            self._dependencies[depname] = []
            self._full_dependencies[depname] = []
        #self._reverse_dependencies = {} # TODO: Need this?
        
        self._plotter = Plotter()
        self._saver = Saver()
        self._planner = Planner()
        
        
        """
        # Plan of what to compute now and in near future
        self._plan = defaultdict(lambda: defaultdict(int))
        
        # Keep track of which fields have been finalized
        self._finalized = {}

        # Cache of computed values needed for planned computations
        self._cache = defaultdict(dict)

        # Keep track of how many times .get has called each field.compute, for administration:
        self._compute_counts = defaultdict(int) # Actually only used for triggering "before_first_compute"
        
        # Keep track of how many times update_all has been called
        self._update_all_count = 0

        # Keep track of last (time, timestep) computation of each field was triggered directly
        self._last_trigger_time = defaultdict(lambda: (-1e16,-1e16))

        # Cache for plotting
        self._plot_cache = {}

        # Callback to be called with fields where the 'callback' action is enabled
        # Signature: ppcallback(field, data, t, timestep)
        self._callback = None

        # Hack to make these objects available throughout during update... Move these to a struct?
        self._problem = None
        self._spaces = None
        self._solution = None
        
        self._timer = Timer()
        """

    @classmethod
    def default_params(cls):
        params = ParamDict(
            casedir=".",
            enable_timer=False,
            extrapolate=True,
            )
        return params

    def _insert_in_sorted_list(self, fieldname):
        # Topological ordering of all fields, so that all dependencies are taken care of

        # If already in list, assuming there's no change to dependencies
        if fieldname in self._sorted_fields_keys:
            return

        # Find largest index of dependencies in sorted list
        deps = [dep[0] for dep in self._dependencies[fieldname] if dep[0] not in builtin_fields]
        max_index = max([-1]+[self._sorted_fields_keys.index(dep) for dep in deps])

        # Insert item after all its dependencies
        self._sorted_fields_keys.insert(max_index+1, fieldname)

    def _find_dependencies(self, field):
        "Read dependencies from source code in field.compute function"
        
        # Get source of compute and after_last_compute
        s = inspect.getsource(field.compute)
        s += inspect.getsource(field.after_last_compute)
        s = strip_code(s) # Removes all comments, empty lines etc.

        # Remove comment blocks
        s = s.split("'''")
        s = s[0::2]
        s = ''.join(s)
        s = s.split('"""')
        s = s[0::2]
        s = ''.join(s)
        s = strip_code(s)

        # Get argument names for the compute function
        args = inspect.getargspec(field.compute)[0]
        self_arg = args[0]
        pp_arg = args[1]

        # Read the code for dependencies
        deps = []
        deps_raw = re.findall(pp_arg+".get\((.+)\)", s)
        for dep in deps_raw:
            # Split into arguments (name, timestep)
            dep = dep.split(',')

            # Append default 0 if dependent timestep not specified
            if len(dep) == 1:
                dep.append(0)

            # Convert timestep to int
            dep[1] = int(dep[1])

            # Get field name from string literal or through string variable
            dep[0] = dep[0].strip(' ').replace('"', "'")
            if "'" in dep[0]:
                # If pp.get('Velocity')
                dep[0] = dep[0].replace("'","")
            else:
                # If pp.get(self.somevariable), get the string hiding at self.somevariable
                dep[0] = eval(dep[0].replace(self_arg, "field", 1))

                # TODO: Test alternative solution without eval (a matter of principle) and with better check:
                #s, n = dep[0].split(".")
                #assert s == self_arg, "Only support accessing strings through self."
                #dep[0] = getattr(field, n)

                # TODO: Possible to get variable in other cases through more introspection?
                #       Probably not necessary, just curious.

            # Append to dependencies
            deps.append(tuple(dep))
            
        # Make unique (can happen that deps are repeated in rare cases)
        return sorted(set(deps))

    def add_field(self, field):
        "Add field to postprocessor. Recursively adds basic dependencies."
        # Did we get a field name instead of a field?
        if isinstance(field, str):
            if field in builtin_fields:
                return None
            elif field in self._fields.keys():
                # Field of this name already exists, no need to add it again
                return self._fields[field]
            elif field in basic_fields:
                # Create a proper field object from known field name with negative end time,
                # so that it is never triggered directly
                field = field_classes[field](params={"end_time":-1e16, "end_timestep": -1e16}) 
            elif field in meta_fields:
                error("Meta field %s cannot be constructed by name because the field instance requires parameters." % field)
            else:
                error("Unknown field name %s" % field)

        # Note: If field already exists, replace anyway to overwrite params, this
        # typically happens when a fields has been created implicitly by dependencies.
        # This is a bit unsafe though, the user might add a field twice with different parameters...
        # Check that at least the same name is not used for different field classes:
        assert type(field) == type(self._fields.get(field.name,field))
        
        # Add fields explicitly specified by field
        self.add_fields(field.add_fields())

        # Analyze dependencies of field through source inspection
        deps = self.find_dependencies(field)

        # Add dependent fields to self._fields (this will add known fields by name)
        for depname in set(d[0] for d in deps) - set(self._fields.keys()):
            self.add_field(depname)

        # Build full dependency list
        full_deps = []
        existing_full_deps = set()
        for dep in deps:
            depname, ts = dep
            for fdep in self._full_dependencies[depname]:
                # Sort out correct (relative) timestep of dependency
                fdepname, fts = fdep
                fts += ts
                fdep = (fdepname, fts)
                if fdep not in existing_full_deps:
                    existing_full_deps.add(fdep)
                    full_deps.append(fdep)
            existing_full_deps.add(dep)
            full_deps.append(dep)
        
        # Add field to internal data structures
        self._fields[field.name] = field
        self._dependencies[field.name] = deps
        self._full_dependencies[field.name] = full_deps
        self._insert_in_sorted_list(field.name)

        # Returning the field object is useful for testing
        return field

    def add_fields(self, fields):
        "Add several fields at once."
        return [self.add_field(field) for field in fields]
    
    def get(self, name, timestep=0, compute=True, finalize=False):
        """Get the value of a named field at a particular timestep.

        The timestep is relative to now.
        Values are computed at first request and cached.
        """
        cbcflow_log(20, "Getting: %s, %d (compute=%s, finalize=%s)" %(name, timestep, compute, finalize))
        
        # Check cache
        c = self._cache[timestep]
        data = c.get(name, "N/A")
        
        # Check if field has been finalized, and if so,
        # return finalized value
        if name in self._finalized and data == "N/A":
            if compute:
                cbcflow_warning("Field %s has already been finalized. Will not call compute on field." %name)
            return self._finalized[name]
        
        # Hack to access the spaces and problem arguments to update()
        spaces = self._spaces
        problem = self._problem
        
        # Are we attempting to get value from before update was started?
        # Use constant extrapolation if allowed.
        if abs(timestep) > self._update_all_count and data == "N/A":
            if self.params.extrapolate:
                cbcflow_log(20, "Extrapolating %s from %d to %d" %(name, timestep, -self._update_all_count))
                data = self.get(name, -self._update_all_count, compute, finalize)
                c[name] = data
            else:
                raise RuntimeError("Unable to get data from before update was started. \
                                   (%s, timestep: %d, update_all_count: %d)" %(name, timestep, self._update_all_count))

        # Cache miss?
        if data == "N/A":
            if timestep == 0:
                # Ensure before_first_compute is always called once initially
                field = self._fields[name]
                if self._compute_counts[field.name] == 0:
                    init_data = field.before_first_compute(self, spaces, problem)
                    if init_data is not None and field.params["save"]:
                        self._init_metadata_file(name, init_data)
                    self._timer.completed("PP: before first compute %s" %name)

                # Compute value
                if name in self._solution:
                    data = field.convert(self, spaces, problem)
                    self._timer.completed("PP: convert %s" %name)
                else:
                    if compute:
                        data = field.compute(self, spaces, problem)
                        self._timer.completed("PP: compute %s" %name)
                    if finalize:
                        finalized_data = field.after_last_compute(self, spaces, problem)
                        if finalized_data not in [None, "N/A"]:
                            data = finalized_data
                        self._finalized[name] = data
                        self._timer.completed("PP: finalize %s" %name)
                self._compute_counts[field.name] += 1

                # Copy functions to avoid storing references to the same function objects at each timestep
                # NB! In other cases we assume that the fields return a new object for every compute!
                # Check first if we actually will cache this object by looking at 'time to keep' in the plan
                if self._plan[0][name] > 0:
                    if isinstance(data, Function):
                        # TODO: Use function pooling to avoid costly allocations?
                        data = Function(data)

                # Cache it!
                c[name] = data
            else:
                # Cannot compute missing value from previous timestep,
                # dependency handling must have failed
                raise DependencyException(name, timestep=timestep)

        return data

    def _execute_plan(self, t, timestep):
        "Check plan and compute fields in plan."
        
        # Initialize cache for current timestep
        assert not self._cache[0], "Not expecting cached computations at this timestep before plan execution!"
        self._cache[0] = {
            "t": t,
            "timestep": timestep,
            }

        # Loop over all planned field computations
        for name in self._sorted_fields_keys:
            if name in self._plan[0]:
                compute = True
            else:
                compute = False
            
            field = self._fields[name]    
            if self._should_finalize_at_this_time(field, t, timestep):
                finalize = True
            else:
                finalize = False
            
            # If neither finalize or compute triggers, continue
            if not (finalize or compute):
                continue
            
            # Execute computation through get call
            data = self.get(name, compute=compute, finalize=finalize)
            
            # Apply action if it was triggered directly this timestep (not just indirectly)
            #if (data is not None) and (self._last_trigger_time[name][1] == timestep):
            if self._last_trigger_time[name][1] == timestep or finalize:
                for action in ["save", "plot", "callback"]:
                    if field.params[action]:
                        self._apply_action(action, field, data)
    
    def _update_cache(self):
        "Update cache, remove what can be removed"
        new_cache = defaultdict(dict)
        # Loop over cache plans for each timestep
        for ts, plan in self._plan.iteritems():
            # Skip whats not computed yet
            if ts > 0:
                continue
            # Only keep what we have planned to cache
            for name, ttk in plan.iteritems():
                if ttk > 0:
                    # Cache should contain old cached values at ts<0 and newly computed values at ts=0
                    data = self._cache[ts].get(name, "N/A")
                    assert data is not "N/A", "Missing cache data!"
                    # Insert data in new cache at position ts-1
                    new_cache[ts-1][name] = data
        self._cache = new_cache


    def update_all(self, solution, t, timestep, spaces, problem):
        "Updates cache, plan, play log and executes plan."

        # TODO: Better design solution to making these variables accessible the right places?
        self._problem = problem
        self._spaces = spaces
        self._solution = solution

        # Update play log
        self._update_play_log(t, timestep)

        # Update cache to keep what's needed later according to plan, forget what we don't need
        self._update_cache()

        # Plan what we need to compute now and in near future based on action triggers and dependencies
        self._update_plan(t, timestep)
        self._timer.completed("PP: updated plan.")

        # Compute what's needed according to plan
        self._execute_plan(t, timestep)
        
        self._update_all_count += 1


    def finalize_all(self, spaces, problem):
        "Finalize all Fields after last timestep has been computed."
        for name in self._sorted_fields_keys:
            field = self._fields[name]
            if field.params.finalize and name not in self._finalized:
                self.get(name, compute=False, finalize=True)
                
            #finalize_data = field.after_last_compute(self, spaces, problem)
            #if finalize_data is not None and field.params["save"]:
            #    self._finalize_metadata_file(name, finalize_data)





