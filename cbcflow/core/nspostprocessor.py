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

from cbcflow.core.paramdict import ParamDict
from cbcflow.core.parameterized import Parameterized
from cbcflow.utils.core.strip_code import strip_code
from cbcflow.utils.common import cbcflow_warning, hdf5_link, safe_mkdir, timeit, on_master_process, in_serial

from cbcflow.fields import field_classes, basic_fields, meta_fields, PPField

from dolfin import Function, MPI, plot, File, HDF5File, XDMFFile, error

import os, re, inspect, pickle, shelve
from collections import defaultdict
from hashlib import sha1
from shutil import rmtree

# TODO: Extract a Plotter class and a Storage class to separate this logic


def disable_plotting():
    "Disable all plotting if we run in parallell."
    if disable_plotting.value == "init":
        if in_serial():
            disable_plotting.value = False
        else:
            cbcflow_warning("Unable to plot dolfin plots in paralell. Disabling.")
            disable_plotting.value = True
    return disable_plotting.value
disable_plotting.value = "init"

def import_pylab():
    "Set up pylab if available."
    if import_pylab.value == "init":
        if disable_plotting():
            import_pylab.value = None
        else:
            try:
                import pylab
                pylab.ion()
                import_pylab.value = pylab
            except:
                cbcflow_warning("Unable to load pylab. Disabling pylab plotting.")
                import_pylab.value = None
    return import_pylab.value
import_pylab.value = "init"

def dolfin_plotting():
    "Enable dolfin plotting if environment variable DISPLAY is set."
    if dolfin_plotting.value == "init":
        if disable_plotting():
            dolfin_plotting.value = False
        else:
            if 'DISPLAY' in os.environ:
                dolfin_plotting.value = True
            else:
                cbcflow_warning("Did not find display. Disabling dolfin plotting.")
                dolfin_plotting.value = False
    return dolfin_plotting.value
dolfin_plotting.value = "init"


class DependencyException(Exception):
    def __init__(self, fieldname=None, dependency=None, timestep=None, original_exception_msg=None):
        message = []
        if fieldname:
            message += ["Dependency/dependencies not found for field %s." % fieldname]
        if dependency:
            message += ["Dependency %s not functioning." % dependency]
        if timestep:
            message += ["Relative timestep is %d. Are you trying to calculate time-derivatives at t=0?" % (fieldname, timestep)]
        if original_exception_msg:
            message += ["\nOriginal exception was: " + original_exception_msg]
        message = ' '.join(message)
        Exception.__init__(self, message)

# Fields available through pp.get(name) even though they have no PPField class
builtin_fields = ("t", "timestep")

class NSPostProcessor(Parameterized):
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

        # Plan of what to compute now and in near future
        self._plan = defaultdict(lambda: defaultdict(int))

        # Cache of computed values needed for planned computations
        self._cache = defaultdict(dict)

        # Keep track of how many times .get has called each field.compute, for administration:
        self._compute_counts = defaultdict(int) # Actually only used for triggering "before_first_compute"

        # Keep track of last (time, timestep) computation of each field was triggered directly
        self._last_trigger_time = defaultdict(lambda: (-1e16,-1e16))

        # Caches for file storage
        self._datafile_cache = {}

        # Cache for plotting
        self._plot_cache = {}

        # Callback to be called with fields where the 'callback' action is enabled
        # Signature: ppcallback(field, data, t, timestep)
        self._callback = None

        # Hack to make these objects available throughout during update... Move these to a struct?
        self._problem = None
        self._spaces = None
        self._solution = None

    @classmethod
    def default_params(cls):
        params = ParamDict(
            casedir=".",
            enable_timer=False,
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

    def find_dependencies(self, field):
        "Read dependencies from source code in field.compute function"
        
        # Get source
        s = inspect.getsource(field.compute)
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
                # Create a proper field object from known field name with default params
                #field = field_classes[field](params=None)
                field = field_classes[field](params={"end_time":-1e16, "end_timestep": -1e16}) # ?Why not None?
            elif field in meta_fields:
                error("Meta field %s cannot be constructed by name because the field instance requires parameters." % field)
            else:
                error("Unknown field name %s" % field)

        # Note: If field already exists, replace anyway to overwrite params, this
        # typically happens when a fields has been created implicitly by dependencies.
        # This is a bit unsafe though, the user might add a field twice with different parameters...
        # Check that at least the same name is not used for different field classes:
        assert type(field) == type(self._fields.get(field.name,field))

        # Analyse dependencies of field through magic
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

    def _should_compute_at_this_time(self, field, t, timestep):
        "Check if field is to be computed at current time"
        # If we later wish to move configuration of field compute frequencies to NSPostProcessor,
        # it's easy to swap here with e.g. fp = self._field_params[field.name]
        fp = field.params

        # Limit by timestep interval
        s = fp.start_timestep
        e = fp.end_timestep
        #if s > timestep or timestep > e:
        if not (s <= timestep <= e):
            return False

        # Limit by time interval
        s = fp.start_time
        e = fp.end_time
        #if s > t or t > e:
        eps = 1e-10
        if not (s-eps <= t <= e+eps):
            return False

        # Limit by frequency (accept if no previous data)
        pct, pcts = self._last_trigger_time[field.name]
        if timestep - pcts < fp.stride_timestep:
            return False
        if t - pct < fp.stride_time:
            return False

        # Accept!
        return True

    def get(self, name, timestep=0):
        """Get the value of a named field at a particular.

        The timestep is relative to now.
        Values are computed at first request and cached.
        """
        # Hack to access the spaces and problem arguments to update()
        spaces = self._spaces
        problem = self._problem

        # Check cache
        c = self._cache[timestep]
        data = c.get(name)

        # Cache miss?
        if data is None:
            if timestep == 0:
                # Ensure before_first_compute is always called once initially
                field = self._fields[name]
                if self._compute_counts[field.name] == 0:
                    init_data = field.before_first_compute(self, spaces, problem)
                    if init_data is not None and field.params["save"]:
                        self._init_metadata_file(name, init_data)

                # Compute value
                t0 = timeit()
                if not name in self._solution:
                    data = field.compute(self, spaces, problem)
                elif self._solution[name]:
                    data = field.convert(self, spaces, problem)
                self._compute_counts[field.name] += 1
                if self.params.enable_timer:
                    timeit(t0, "time to compute %s" % (name,))

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
                raise DependencyException(name, timestep)

        return data

    def _apply_action(self, action, field, data):
        "Apply some named action to computed field data."
        if not field.name in self._cache[0]:
            error("Field '%s' is not in cache, this should not be possible." % field.name)

        t0 = timeit()
        if action == "save":
            self._action_save(field, data)

        elif action == "plot":
            self._action_plot(field, data)

        elif action == "callback":
            self._action_callback(field, data)

        else:
            error("Unknown action %s." % action)
        if self.params.enable_timer:
            timeit(t0, "time to %s %s" % (action, field.name))

    def _action_callback(self, field, data):
        "Apply the 'callback' action to computed field data."
        if callable(self._callback):
            self._callback(field, data, self.get("t"), self.get("timestep"))

    def _get_save_formats(self, field, data):
        if field.params.save_as == PPField.default_save_as():
            # Determine proper file formats from data type if not specifically provided
            if isinstance(data, Function):
                save_as = ['xdmf', 'hdf5']
            elif isinstance(data, (float, int, list, tuple, dict)):
                save_as = ['txt', 'shelve']
            else:
                error("Unknown data type %s, cannot determine file type automatically." % type(data).__name__)
        else:
            if isinstance(field.params.save_as, (list, tuple)):
                save_as = list(field.params.save_as)
            else:
                save_as = [field.params.save_as]
        return save_as

    def get_casedir(self):
        return self.params.casedir

    def _clean_casedir(self):
        "Cleans out all files produced by cbcflow in the current casedir."
        if on_master_process():
            if os.path.isdir(self.get_casedir()):
                playlogfilename = os.path.join(self.get_casedir(), "play.db")
                if os.path.isfile(playlogfilename):
                    playlog = shelve.open(playlogfilename, 'r')

                    all_fields = []
                    for k,v in playlog.items():
                        all_fields += v.get("fields", {}).keys()
    
                    all_fields = list(set(all_fields))
                    playlog.close()
                    
                    for field in all_fields:
                        rmtree(os.path.join(self.get_casedir(), field))
                    
                    for f in ["mesh.hdf5", "play.db", "params.txt", "params.pickle"]:
                        if os.path.isfile(os.path.join(self.get_casedir(), f)):
                            os.remove(os.path.join(self.get_casedir(), f))


    def _create_casedir(self):
        casedir = self.params.casedir
        safe_mkdir(casedir)
        return casedir

    def get_savedir(self, field_name):
        "Returns savedir for given fieldname"
        return os.path.join(self.params.casedir, field_name)

    def _create_savedir(self, field_name):
        self._create_casedir()
        savedir = self.get_savedir(field_name)
        safe_mkdir(savedir)
        return savedir

    def _init_metadata_file(self, field_name, init_data):
        savedir = self._create_savedir(field_name)
        if on_master_process():
            metadata_filename = os.path.join(savedir, 'metadata.db')
            metadata_file = shelve.open(metadata_filename)
            metadata_file["init_data"] = init_data
            metadata_file.close()

    def _finalize_metadata_file(self, field_name, finalize_data):
        if on_master_process():
            savedir = self.get_savedir(field_name)
            metadata_filename = os.path.join(savedir, 'metadata.db')
            metadata_file = shelve.open(metadata_filename)
            metadata_file["finalize_data"] = finalize_data
            metadata_file.close()

    def _update_metadata_file(self, field_name, data, t, timestep, save_as, metadata):
        if on_master_process():
            savedir = self.get_savedir(field_name)
            metadata_filename = os.path.join(savedir, 'metadata.db')
            metadata_file = shelve.open(metadata_filename)

            # Store some data the first time
            if "type" not in metadata_file:
                # Data about type and formats
                metadata_file["type"] = type(data).__name__
                metadata_file["saveformats"] = list(set(save_as+metadata_file.get("saveformats", [])))
                # Data about function space
                if isinstance(data, Function):
                    metadata_file["element"] = repr(data.element(),)
                    metadata_file["element_degree"] = repr(data.element().degree(),)
                    metadata_file["element_family"] = repr(data.element().family(),)
                    metadata_file["element_value_shape"] = repr(data.element().value_shape(),)
            # Store some data each timestep
            metadata_file[str(timestep)] = metadata
            metadata_file[str(timestep)]["t"] = t

            # Flush file between timesteps
            metadata_file.close()

    def _get_datafile_name(self, field_name, saveformat, timestep):
        # These formats produce a new file each time
        counted_formats = ('xml', 'xml.gz')

        metadata = {}

        # Make filename, with or without save count in name
        if saveformat in counted_formats:
            filename = "%s%d.%s" % (field_name, timestep, saveformat)
            # If we have a new filename each time, store the name in metadata
            #metadata = [('filename', filename)]
            metadata['filename'] = filename
        elif saveformat == "shelve":
            filename = "%s.%s" % (field_name, "db")
        else:
            filename = "%s.%s" % (field_name, saveformat)
            if saveformat == 'hdf5':
                metadata['dataset'] = field_name+str(timestep)

        savedir = self.get_savedir(field_name)
        fullname = os.path.join(savedir, filename)
        return fullname, metadata

    def _update_pvd_file(self, field_name, saveformat, data, timestep, t):
        assert isinstance(data, Function)
        assert saveformat == "pvd"
        fullname, metadata = self._get_datafile_name(field_name, saveformat, timestep)
        key = (field_name, saveformat)
        datafile = self._datafile_cache.get(key)
        if datafile is None:
            datafile = File(fullname)
            self._datafile_cache[key] = datafile
        datafile << data
        return metadata

    def _update_xdmf_file(self, field_name, saveformat, data, timestep, t):
        assert isinstance(data, Function)
        assert saveformat == "xdmf"
        fullname, metadata = self._get_datafile_name(field_name, saveformat, timestep)
        key = (field_name, saveformat)
        datafile = self._datafile_cache.get(key)
        if datafile is None:
            datafile = XDMFFile(fullname)
            datafile.parameters["rewrite_function_mesh"] = False
            datafile.parameters["flush_output"] = True
            self._datafile_cache[key] = datafile
        datafile << (data, t)
        return metadata

    def _update_hdf5_file(self, field_name, saveformat, data, timestep, t):
        assert saveformat == "hdf5"
        fullname, metadata = self._get_datafile_name(field_name, saveformat, timestep)
        
        # Create "good enough" hash. This is done to avoid data corruption when restarted from
        # different number of processes, different distribution or different function space
        local_hash= sha1()
        local_hash.update(str(data.function_space().mesh().num_cells()))
        local_hash.update(str(data.function_space().ufl_element()))
        local_hash.update(str(data.function_space().dim()))
        local_hash.update(str(MPI.num_processes()))
        
        # Global hash (same on all processes), 10 digits long
        hash = str(int(MPI.sum(int(local_hash.hexdigest(), 16))%1e10)).zfill(10)
        
        # Open HDF5File
        if not os.path.isfile(fullname):
            datafile = HDF5File(fullname, 'w')
        else:
            datafile = HDF5File(fullname, 'a')
        
        # Write to hash-dataset if not yet done
        if not datafile.has_dataset(hash) or not datafile.has_dataset(hash+"/"+field_name):
            datafile.write(data, str(hash)+"/"+field_name)
            
        if not datafile.has_dataset("Mesh"):
            datafile.write(data.function_space().mesh(), "Mesh")
        
        # Write vector to file
        # TODO: Link vector when function has been written to hash
        datafile.write(data.vector(), field_name+str(timestep)+"/vector")
        del datafile

        # Link information about function space from hash-dataset
        hdf5_link(fullname, str(hash)+"/"+field_name+"/x_cell_dofs", field_name+str(timestep)+"/x_cell_dofs")
        hdf5_link(fullname, str(hash)+"/"+field_name+"/cell_dofs", field_name+str(timestep)+"/cell_dofs")
        hdf5_link(fullname, str(hash)+"/"+field_name+"/cells", field_name+str(timestep)+"/cells")

        return metadata

    def _update_xml_file(self, field_name, saveformat, data, timestep, t):
        assert saveformat == "xml"
        fullname, metadata = self._get_datafile_name(field_name, saveformat, timestep)
        datafile = File(fullname)
        datafile << data
        return metadata

    def _update_xml_gz_file(self, field_name, saveformat, data, timestep, t):
        assert saveformat == "xml.gz"
        fullname, metadata = self._get_datafile_name(field_name, saveformat, timestep)
        datafile = File(fullname)
        datafile << data
        return metadata

    def _update_txt_file(self, field_name, saveformat, data, timestep, t):
        # TODO: Identify which more well defined data formats we need
        assert saveformat == "txt"
        fullname, metadata = self._get_datafile_name(field_name, saveformat, timestep)
        if on_master_process():
            datafile = open(fullname, 'a')
            datafile.write(str(data))
            datafile.write("\n")
            datafile.close()
        return metadata

    def _update_shelve_file(self, field_name, saveformat, data, timestep, t):
        assert saveformat == "shelve"
        fullname, metadata = self._get_datafile_name(field_name, saveformat, timestep)
        if on_master_process():
            datafile = shelve.open(fullname)
            datafile[str(timestep)] = data
            datafile.close()
            
        return metadata

    def _fetch_play_log(self):
        casedir = self.get_casedir()
        play_log_file = os.path.join(casedir, "play.db")
        play_log = shelve.open(play_log_file)
        return play_log

    def _update_play_log(self, t, timestep):
        if on_master_process():
            play_log = self._fetch_play_log()
            if str(timestep) in play_log:
                play_log.close()
                return
            play_log[str(timestep)] = {"t":float(t)}
            play_log.close()

    def _fill_play_log(self, field, timestep, save_as):
        if on_master_process():
            play_log = self._fetch_play_log()
            timestep_dict = dict(play_log[str(timestep)])
            if "fields" not in timestep_dict:
                timestep_dict["fields"] = {}
            timestep_dict["fields"][field.name] = {"type": field.__class__.shortname(), "save_as": save_as}
            play_log[str(timestep)] = timestep_dict
            play_log.close()

    def store_params(self, params):
        "Store parameters in casedir as params.pickle and params.txt."
        casedir = self._create_casedir()

        pfn = os.path.join(casedir, "params.pickle")
        with open(pfn, 'w') as f:
            pickle.dump(params, f)

        tfn = os.path.join(casedir, "params.txt")
        with open(tfn, 'w') as f:
            f.write(str(params))

    def store_mesh(self, mesh):
        "Store mesh in casedir to mesh.hdf5 (dataset Mesh) in casedir."
        casedir = self.get_casedir()
        meshfile = HDF5File(os.path.join(casedir, "mesh.hdf5"), 'w')
        meshfile.write(mesh, "Mesh")
        del meshfile

    def _action_save(self, field, data):
        "Apply the 'save' action to computed field data."
        field_name = field.name

        # Create save folder first time
        self._create_savedir(field_name)

        # Get current time (assuming the cache contains
        # valid 't' and 'timestep' at each step)
        t = self.get("t")
        timestep = self.get('timestep')

        # Get list of file formats
        save_as = self._get_save_formats(field, data)

        # Collect metadata shared between data types
        metadata = {
            'timestep': timestep,
            'time': t,
            }

        # Rename Functions to get the right name in file
        # (NB! This has the obvious side effect!)
        # TODO: We don't need to cache a distinct Function
        # object like we do for plotting, or?
        if isinstance(data, Function):
            data.rename(field_name, "Function produced by cbcflow postprocessing.")

        # Write data to file for each filetype
        for saveformat in save_as:
            # Write data to file depending on type
            if saveformat == 'pvd':
                metadata[saveformat] = self._update_pvd_file(field_name, saveformat, data, timestep, t)
            elif saveformat == 'xdmf':
                metadata[saveformat] = self._update_xdmf_file(field_name, saveformat, data, timestep, t)
            elif saveformat == 'xml':
                metadata[saveformat] = self._update_xml_file(field_name, saveformat, data, timestep, t)
            elif saveformat == 'xml.gz':
                metadata[saveformat] = self._update_xml_gz_file(field_name, saveformat, data, timestep, t)
            elif saveformat == 'txt':
                metadata[saveformat] = self._update_txt_file(field_name, saveformat, data, timestep, t)
            elif saveformat == 'hdf5':
                metadata[saveformat] = self._update_hdf5_file(field_name, saveformat, data, timestep, t)
            elif saveformat == 'shelve':
                metadata[saveformat] = self._update_shelve_file(field_name, saveformat, data, timestep, t)
            else:
                error("Unknown save format %s." % (saveformat,))

        # Write new data to metadata file
        self._update_metadata_file(field_name, data, t, timestep, save_as, metadata)
        
        self._fill_play_log(field, timestep, save_as)


    def _action_plot(self, field, data):
        "Apply the 'plot' action to computed field data."
        if disable_plotting():
            return
        if isinstance(data, Function):
            self._plot_dolfin(field.name, data)
        elif isinstance(data, float):
            self._plot_pylab(field.name, data)
        else:
            cbcflow_warning("Unable to plot object %s of type %s." % (field.name, type(data)))

    def _plot_dolfin(self, field_name, data):
        "Plot field using dolfin plot command"
        if not dolfin_plotting():
            return

        # Get current time
        t = self.get("t")
        timestep = self.get('timestep')

        # Plot or re-plot
        plot_object = self._plot_cache.get(field_name)
        if plot_object is None:
            plot_object = plot(data, title=field_name, **self._fields[field_name].params.plot_args)
            self._plot_cache[field_name] = plot_object
        else:
            plot_object.plot(data)

        # Set title and show
        title = "%s, t=%0.4g, timestep=%d" % (field_name, t, timestep)
        plot_object.parameters["title"] = title

    def _plot_pylab(self, field_name, data):
        "Plot using pylab if field is a single scalar."
        pylab = import_pylab()
        if not pylab:
            return

        # Hack to access the spaces and problem arguments to update()
        problem = self._problem

        # Get current time
        t = self.get("t")
        timestep = self.get('timestep')

        # Values to plot
        x = t
        y = data

        # Plot or re-plot
        plot_data = self._plot_cache.get(field_name)
        if plot_data is None:
            figure_number = len(self._plot_cache)
            pylab.figure(figure_number)

            xdata = [x]
            ydata = [y]
            newmin = min(ydata)
            newmax = max(ydata)

            plot_object, = pylab.plot(xdata, ydata)
            self._plot_cache[field_name] = plot_object, figure_number, newmin, newmax
        else:
            plot_object, figure_number, oldmin, oldmax = plot_data
            pylab.figure(figure_number)

            xdata = list(plot_object.get_xdata())
            ydata = list(plot_object.get_ydata())
            xdata.append(x)
            ydata.append(y)
            newmin = min(ydata)
            newmax = max(ydata)

            # Heuristics to avoid changing axis bit by bit, which results in fluttering plots
            # (Based on gut feeling, feel free to adjust these if you have a use case it doesnt work for)
            if newmin < oldmin:
                # If it has decreased, decrease by at least this factor
                #ymin = min(newmin, oldmin*0.8) # TODO: Negative numbers?
                ymin = newmin
            else:
                ymin = newmin
            if newmax > oldmax:
                # If it has increased, increase by at least this factor
                #ymax = max(newmax, oldmax*1.2) # TODO: Negative numbers?
                ymax = newmax
            else:
                ymax = newmax

            # Need to store min/max for the heuristics to work
            self._plot_cache[field_name] = plot_object, figure_number, ymin, ymax

            plot_object.set_xdata(xdata)
            plot_object.set_ydata(ydata)

            pylab.axis([problem.params.T0, problem.params.T, ymin, ymax])

        # Set title and show
        title = "%s, t=%0.4g, timestep=%d, min=%.2g, max=%.2g" % (field_name, t, timestep, newmin, newmax)
        plot_object.get_axes().set_title(title)
        pylab.xlabel("t")
        pylab.ylabel(field_name)
        pylab.draw()

    def _rollback_plan(self, t, timestep):
        "Roll plan one timestep and countdown how long to keep stuff."
        tss = sorted(self._plan.keys())
        new_plan = defaultdict(lambda: defaultdict(int))

        # Countdown time to keep each data item and only keep what we still need
        for ts in tss:
            for name, ttk in self._plan[ts].items():
                if ttk > 0:
                    new_plan[ts-1][name] = ttk - 1
        self._plan = new_plan

    def _update_plan(self, t, timestep):
        "Update plan for new timestep."
        # ttk = timesteps to keep
        #self._plan[-1][name] = ttk # How long to cache what's already computed
        #self._plan[0][name] = ttk  # What to compute now and how long to cache it
        #self._plan[1][name] = ttk  # What to compute in future and how long to cache it

        self._rollback_plan(t, timestep)

        # Loop over all fields that are triggered for computation at this timestep
        triggered_fields = [(name, field) for name, field in self._fields.iteritems()
                            if self._should_compute_at_this_time(field, t, timestep)]

        for name, field in triggered_fields:
            deps = self._full_dependencies[name]
            if deps:
                # Need to plan ahead this many steps (ts is non-positive)
                offset = abs(min(ts for depname, ts in deps))
                # Plan computation of dependenices:
                for depname, ts in deps:
                    # Store how long we need to cache this computation, highest of offset and old plan
                    oldttk = self._plan[ts+offset].get(depname, 0)
                    ttk = max(oldttk, offset)
                    self._plan[ts+offset][depname] = ttk
            else:
                # No planning ahead, compute right away
                offset = 0

            # Plan computation of this field:
            oldttk = self._plan[offset].get(name, 0)
            ttk = max(oldttk, offset)
            self._plan[offset][name] = ttk

            # Store compute trigger times to keep track of compute intervals
            self._last_trigger_time[field.name] = (t, timestep)

    def _execute_plan(self, t, timestep):
        "Check plan and compute fields in plan."
        # Initialize cache for current timestep
        assert not self._cache[0], "Not expecting cached computations at this timestep before plan execution!"
        self._cache[0] = {
            "t": t,
            "timestep": timestep,
            }
        # Loop over all planned field computations
        fields_to_compute = [name for name in self._sorted_fields_keys if name in self._plan[0]]
        for name in fields_to_compute:
            field = self._fields[name]
            # Execute computation through get call
            try:
                data = self.get(name)
            except DependencyException as e:
                cbcflow_warning(e.message)
                data = None

            # Apply action if it was triggered directly this timestep (not just indirectly)
            if (data is not None) and (self._last_trigger_time[name][1] == timestep):
                for action in ["save", "plot", "callback"]:
                    if field.params[action]:
                        self._apply_action(action, field, data)

        if 0: # Debugging:
            s1 = set(fields_to_compute)
            s2 = set(self._cache[0])
            s2.remove("t")
            s2.remove("timestep")
            s2m1 = s2 - s1
            s1m2 = s1 - s2
            if s2m1 or s1m2:
                print "\nCompute lists differ:"
                print sorted(s2m1) #s2-s1 = fields in s2 but not in s1, that is in expected cache but not in computation list
                print sorted(s1m2)
                print "Plan is:"
                print self._plan[0]
                print "Sorted field keys:"
                print self._sorted_fields_keys
                print "fields_to_compute:"
                print fields_to_compute
                print "cache[0]:"
                print self._cache[0]
                print
        
        # Wrong: Try adding e.g. TimeDerivative("Stress") only to a postprocessor. This check fails, but planning algorithm seems correct
        #assert len(fields_to_compute) == len(self._cache[0])-2, "This should hold if planning algorithm works properly?"
        

    def _update_cache(self):
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
                    data = self._cache[ts].get(name)
                    assert data is not None, "Missing cache data!"
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

        # Compute what's needed according to plan
        self._execute_plan(t, timestep)

        # Reset hack to make these objects available throughout during update, to
        # make sure these objects are not referenced after this update call is over
        #self._problem = None
        #self._spaces = None
        #self._solution = None

        # FIXME: Tests are calling pp.get after update_all, this fails for two reasons:
        # - the hack reset above
        # - the cache rotation above (previously cache rotation was at the beginning of update_all)
        # We can't just use pp.get(name,-1) from the tests either, because that will not
        # allow new computations. Do we want to allow stand-alone get? Find a solution!
        # Probably involves updating cache at beginning, so we don't forget data until the
        # top of the next update_all call.

    def finalize_all(self, spaces, problem):
        "Finalize all PPFields after last timestep has been computed."
        for name in self._sorted_fields_keys:
            field = self._fields[name]
            finalize_data = field.after_last_compute(self, spaces, problem)
            if finalize_data is not None and field.params["save"]:
                self._finalize_metadata_file(name, finalize_data)
