from .PPField import PPField
from ..utils import *
from dolfin import *

class SubFunction(PPField):
    def __init__(self, field, submesh, params=None, label=None):
        PPField.__init__(self, params, label)
        
        import imp
        try:
            imp.find_module("mpi4py")
        except:
            raise ImportError("Can't find module mpi4py. This is required for SubFunction.")

        self.submesh = submesh
        
        # Store only name, don't need the field
        if isinstance(field, PPField):
            value = field.name
        self.valuename = field

    @property
    def name(self):
        n = "SubFunction_%s" % self.valuename
        if self.label: n += "_"+self.label
        return n
    
    
    def before_first_compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)
        
        V = u.function_space()
        element = V.ufl_element()        
        family = element.family()
        degree = element.degree()
        
        if u.rank() == 1:
            FS = VectorFunctionSpace(self.submesh, family, degree)
            FS_scalar = FS.sub(0).collapse()
            self.u0 = Function(FS_scalar)
            self.u1 = Function(FS_scalar)
            self.u2 = Function(FS_scalar)
            
        elif u.rank() == 0:
            FS = FunctionSpace(self.submesh, family, degree)
        else:
            raise Exception("Does not support TensorFunctionSpace yet")
        
        self.u = Function(FS, name=self.name)


    def compute(self, pp, spaces, problem):       
        u = pp.get(self.valuename)
        
        if u.rank() == 1:
            u0, u1, u2 = u.split()
            
            self.u0.assign(interpolate_nonmatching_mesh(u0, self.u0))
            self.u1.assign(interpolate_nonmatching_mesh(u1, self.u1))
            self.u2.assign(interpolate_nonmatching_mesh(u2, self.u2))
            
            self.u.assign(project(as_vector([self.u0, self.u1, self.u2])))           
        elif u.rank() == 0:
            self.u.assign(interpolate_nonmatching_mesh(u, self.u))
            
        return self.u
