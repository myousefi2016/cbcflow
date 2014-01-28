
from fractions import gcd
from os.path import isfile, isdir, join
from os import listdir, makedirs

from dolfin import Function, File, MPI, TestFunction, assemble, inner, dx, project, HDF5File
import shelve   
from ...core.paramdict import ParamDict
from ...core.parameterized import Parameterized

from ...core.utils import cbcflow_warning

class PPField(Parameterized):
    def __init__(self, params=None, label=None):
        Parameterized.__init__(self, params)
        if label:
            self.label = str(label)
        else:
            self.label = None

    # --- Parameters

    @classmethod
    def default_save_as(cls):
        return "determined by data"

    @classmethod
    def default_params(cls):
        params = ParamDict(
            # Configure direct compute requests through timestep counting
            start_timestep = -1e16,
            end_timestep = 1e16,
            stride_timestep = 1,

            # Configure direct compute requests through physical time intervals
            start_time = -1e16,
            end_time = 1e16,
            stride_time = 1e-16,

            # Trigger action after each direct compute request
            plot = False,
            save = False,
            callback = False,

            # Configure computing
            project = False, # This is the safest approach
            assemble = True, # This is faster but only works for for DG0
            interpolate = False, # This will be the best when properly implemented in fenics

            # Configure saving
            save_as = cls.default_save_as(),

            # Configure plotting
            plot_args={},
            
            # Solution switch
            is_solution = False,
            )
        return params

    @property
    def name(self):
        "Return name of field, by default the classname but can be overloaded in subclass."
        n = self.__class__.__name__
        if self.label: n += "_"+self.label
        return n

    # --- Main interface

    def before_first_compute(self, pp, spaces, problem):
        "Called prior to the simulation timeloop."
        pass

    def after_last_compute(self, pp, spaces, problem):
        "Called after the simulation timeloop."
        pass

    def compute(self, pp, spaces, problem):
        "Called each time the quantity should be computed."
        raise NotImplementedError("A PPField must implement the compute function!")

    def convert(self, pp, spaces, problem):
        """Called if quantity is input to NSPostProcessor.update_all"""
        
        # Load data from disk (this is used in replay functionality)
        # The structure of the dict pp._solution[self.name] is determined in nsreplay.py
        if isinstance(pp._solution[self.name], dict):
            saveformat = pp._solution[self.name]["format"]
            if saveformat == 'hdf5':
                hdf5filepath = join(pp._get_casedir(), self.name, self.name+".hdf5")
                hdf5file = HDF5File(hdf5filepath, 'r')
                dataset = self.name+str(pp._solution[self.name]["save_count"])
                hdf5file.read(pp._solution[self.name]["function"], dataset)
                pp._solution[self.name] = pp._solution[self.name]["function"]
            elif saveformat in ["xml", "xml.gz"]:
                save_count = pp._solution[self.name]["save_count"]
                xmlfilename = self.name+str(save_count)+"."+saveformat
                xmlfilepath = join(pp._get_casedir(), self.name, xmlfilename)
                function = pp._solution[self.name]["function"]
                function.assign(Function(function.function_space(), xmlfilepath))
                pp._solution[self.name] = pp._solution[self.name]["function"]
            elif saveformat == "shelve":
                shelvefilepath = join(pp._get_casedir(), self.name, self.name+".db")
                shelvefile = shelve.open(shelvefilepath)
                save_count = pp._solution[self.name]["save_count"]
                pp._solution[self.name] = shelvefile[str(save_count)]

        return pp._solution[self.name]

    # --- Helper functions

    def expr2function(self, expr, function):

        space = function.function_space()

        if self.params.assemble:
            # Compute average values of expr for each cell and place in a DG0 space

            # TODO: Get space from pool
            #shape = expr.shape()
            #space = pp.space_pool.get_custom_space("DG", 0, shape)
            #target = pp.function_pool.borrow_function(space)

            test = TestFunction(space)
            scale = 1.0 / space.mesh().ufl_cell().volume
            assemble(scale*inner(expr, test)*dx(), tensor=function.vector())
            return function

        elif self.params.project:
            # TODO: Avoid superfluous function creation by allowing project(expr, function=function) or function.project(expr)
            function.assign(project(expr, space))
            return function

        elif self.params.interpolate:
            # TODO: Need interpolation with code generated from expr, waiting for uflacs work.
            function.interpolate(expr) # Currently only works if expr is a single Function
            return function

        else:
            error("No action selected, need to choose either assemble, project or interpolate.")