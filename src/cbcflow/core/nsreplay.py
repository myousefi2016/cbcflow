__author__ = "Oeyvind Evju <oyvinev@simula.no>"
__date__ = "2013-10-18"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"


from .parameterized import Parameterized
from .paramdict import ParamDict
from .nsproblem import NSProblem
from .nspostprocessor import NSPostProcessor
from ..postprocessing.PPField import PPField
from ..postprocessing import field_classes
from .spaces import NSSpacePoolSplit
import pickle
import os
from dolfin import HDF5File, Mesh, Function, FunctionSpace, VectorFunctionSpace, TensorFunctionSpace, BoundaryMesh
from .utils import cbcflow_print
import inspect, shelve
from .utils import cbcflow_warning


class NSReplay(Parameterized):
    """ Replay class for postprocessing exisiting solution data. """
    def __init__(self, postprocessor, params=None):
        Parameterized.__init__(self, params)
        self.postproc = postprocessor
        self.data_dict = {}
        
        self._functions = {}
        self._hdf5_files = {}
        
    @classmethod
    def default_params(cls):
        params = ParamDict(
            mesh_file = '',
            )
        return params

    def _fetch_history(self):
        data = {}
        #assert(os.path.isfile(os.path.join(self.postproc._get_casedir(), "play.shelve")))
        
        # Read play.shelve
        #play = shelve.open(os.path.join(self.postproc._get_casedir(), "play.shelve"))
        #for key, value in play.items():
        #    data[int(key)] = {"t": value}

        # Read all field metadata
        all_folders = [x[0] for x in os.walk(self.postproc._get_casedir())]
        folders_match = [self.postproc._get_savedir(fieldname) for fieldname in field_classes.keys()]
        
        for d in all_folders:
            if any([f in d for f in folders_match]):
                fieldname = os.path.split(d)[1]
                
                # FIXME: Allow for all fields
                if not fieldname in ["Velocity", "Pressure"]:
                    continue
                metadata_filename = os.path.join(d, "metadata.db")
                if os.path.isfile(metadata_filename):
                    metadata_file = shelve.open(metadata_filename, 'r')
                    for key, value in metadata_file.items():
                        try:
                            key = int(key)
                        except:
                            continue
                        
                        if key not in data:
                            data[key] = {}
                        #if int(key) not in data:
                        #    cbcflow_warning("Unable to find timestep %d in play-file, but found in metadata for field %s." %(int(key), fieldname))
                        #    continue
                        if "hdf5" in value:
                            value['hdf5']['filename'] = fieldname+".hdf5"
                        
                        data[int(key)][fieldname] = value
                        if 't' not in data and 'time' in value:
                            data[int(key)]['t'] = value['time']
        return data
   
       
    def _check_field_coverage(self, plan, fieldname):
        "Find which timesteps fieldname can be computed at"
        timesteps = []
        for ts in plan.keys():
            if self._recursive_dependency_check(plan, ts, fieldname):
                timesteps.append(ts)

        return timesteps
    
        
    def _recursive_dependency_check(self, plan, key, fieldname):
        "Check if field or dependencies exist in plan"        
        if key not in plan:
            return False
        if fieldname == "t":
            return True
        
        fielddata = plan[key].get(fieldname)
        # Return True if data present in plan as a readable data format
        #dataformats = ["hdf5", "xml", "xml.gz", "float", "int", "list", "tuple", "dict"]
        dataformats = ["hdf5"]
        if fielddata:
            if any(d in fielddata for d in dataformats):
                return True
        
        # Check if data should be available in a readable format
        field = self.postproc._fields.get(fieldname)
        if field and field.params.is_solution:
            return False
        
        # If no dependencies, or if not all dependencies exist in plan, return False
        dependencies = self.postproc._dependencies[fieldname]
        if len(dependencies) == 0:
            return False
        else:
            checks = []
            for dep_field, dep_time in dependencies:
                if dep_time != 0:
                    continue
                checks.append(self._recursive_dependency_check(plan, key, dep_field))

            return all(checks)
       
    def _create_function_from_metadata(self, fieldname):
        metadata_filename = os.path.join(self.postproc._get_casedir(), fieldname, "metadata.db")
        if not os.path.isfile(metadata_filename):
            return None
        
        metadata_file = shelve.open(metadata_filename, 'r')

        # Should avoid this, but BoundaryMesh has some bug
        if fieldname == "WSS":
            meshfilename = os.path.join(self.postproc._get_casedir(), 'mesh.hdf5')
            assert os.path.isfile(meshfilename)
            hdf5file = HDF5File(meshfilename, 'r')
            bdry = Mesh()
            hdf5file.read(bdry, "BoundaryMesh")
            del hdf5file

            V = VectorFunctionSpace(bdry, "DG", 0)
            
        else:
            # Create function space from exisiting metadata
            shape = eval(metadata_file["element_value_shape"])
            degree = eval(metadata_file["element_degree"])
            family = eval(metadata_file["element_family"])
            V = self.spaces.get_space(degree, len(shape), family)

        # Create function
        u = Function(V, name=fieldname)
        return u

       
       
    def _get_function(self, fieldname):
        if fieldname not in self._functions:
            self._functions[fieldname] = self._create_function_from_metadata(fieldname)
        return self._functions[fieldname]
       
    def replay(self):
        "Replay problem with given postprocessor."

        # Initiate problem
        paramfile = open(os.path.join(self.postproc._get_casedir(), "params.pickle"), 'rb')
        params = pickle.load(paramfile)
        
        # FIXME: This is ugly, and shouldn't be necessary
        params.problem.num_periods = None
        problem = NSProblem(params.problem)

        # Read mesh
        meshfilename = os.path.join(self.postproc._get_casedir(), "Pressure/Pressure.hdf5")
        assert os.path.isfile(meshfilename), "Unable to find meshfile!"
        meshfile = HDF5File(meshfilename, 'r')
        mesh = Mesh()
        meshfile.read(mesh, "Mesh")
        del meshfile
        
        problem.mesh = mesh
        
        # Set up for replay
        replay_plan = self._fetch_history()
        postprocessors = [] 
        for fieldname, field in self.postproc._fields.items():
            # Check timesteps covered by current field
            keys = self._check_field_coverage(replay_plan, fieldname)
            
            # Get the time dependency for the field
            t_dep = min([dep[1] for dep in self.postproc._dependencies[fieldname]]+[0])
            
            
            # Append field to correct postprocessor
            added_to_postprocessor = False
            for i, (ppkeys, ppt_dep, pp, ppfields) in enumerate(postprocessors):               
                if t_dep == ppt_dep:
                    if t_dep == 0:
                        ppkeys = list(set(keys+ppkeys))
                        ppfields.append(field)
                        added_to_postprocessor = True
                    elif set(keys).issubset(set(ppkeys)):
                        ppfields.append(field)
                        added_to_postprocessor = True

            # Create new postprocessor if no suitable postprocessor found
            if not added_to_postprocessor:
                pp = NSPostProcessor({"casedir": self.postproc._get_casedir()})
                postprocessors.append((keys, t_dep, pp, [field]))

        # Add all fields at once to each postprocessor, to avoid dependency bugs
        for ppkeys, ppt_dep, pp, ppfields in postprocessors:
            pp.add_fields(ppfields)
        
        # Create spaces object
        self.spaces = NSSpacePoolSplit(mesh, params.scheme.u_degree, params.scheme.p_degree)
        
        # Run replay
        for timestep in sorted(replay_plan.keys()):
            #cbcflow_print("Processing timestep %d of %d. %.3f%% complete." %(timestep, len(replay_plan.keys())-1, 100.0*(timestep)/(len(replay_plan.keys())-1)))
            cbcflow_print("Processing timestep %d of %d. %.3f%% complete." %(timestep, max(replay_plan.keys()), 100.0*(timestep)/(max(replay_plan.keys()))))

            # Load solution at this timestep (all available fields)
            saved_fields = replay_plan[timestep]
            solution = {}
            for fieldname in saved_fields:
                if not isinstance(saved_fields[fieldname], dict):
                    continue
                
                if "hdf5" in saved_fields[fieldname]:
                    function = self._get_function(fieldname)
                    dataset = saved_fields[fieldname]["hdf5"]["dataset"]
                    hdf5filename = os.path.join(self.postproc._get_casedir(), fieldname, saved_fields[fieldname]["hdf5"]["filename"])
                    solution[fieldname] = {"HDF5": (HDF5File(hdf5filename, 'r'), dataset, function)}
                #TODO: Fix support for other formats
                """
                elif "xml" in saved_fields[fieldname]:
                    function = self._get_function(fieldname)
                    xmlfilename = os.path.join(self.postproc._get_casedir(), fieldname, saved_fields[fieldname]["xml"]["filename"])
                    solution[fieldname] = {"xml": (xmlfilename, function)}
                elif "xml.gz" in saved_fields[fieldname]:
                    function = self._get_function(fieldname)
                    xmlfilename = os.path.join(self.postproc._get_casedir(), fieldname, saved_fields[fieldname]["xml.gz"]["filename"])
                    solution[fieldname] = {"xml.gz": (xmlfilename, function)}
                elif "shelve" in saved_fields[fieldname]:
                    # Read this data directly
                    datafilename = os.path.join(self.postproc._get_casedir(), fieldname, saved_fields[fieldname]["shelve"]["filename"])
                    datafile = shelve.open(datafilename)
                    solution[fieldname] = datafile[str(timestep)]
                """
            # Cycle through postprocessors and update if required
            for ppkeys, ppt_dep, pp, ppfields in postprocessors:
                if timestep in ppkeys:
                    pp.update_all(solution, saved_fields["t"], timestep, self.spaces, problem)

                    # Update solution to avoid re-reading data
                    solution = pp._solution
                        
                
   

