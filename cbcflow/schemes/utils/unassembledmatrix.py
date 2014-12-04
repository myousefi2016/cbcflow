import os
import inspect
import re
from dolfin import compile_extension_module

def _compile_code():
    def strip_essential_code(filename):
        f = open(filename, 'r').read()
        code = f[f.find("namespace dolfin\n{\n"):f.find("#endif")]
        return code
    
    folder = os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0])
    header = os.path.join(folder, "UnassembledMatrix.h")
    source = "UnassembledMatrix.cpp"
    code = strip_essential_code(header)
    
    try:
        compiled_module = compile_extension_module(code=code, source_directory=os.path.abspath(folder),
                                                sources=[source], include_dirs=[".", os.path.abspath(folder)],
                                                cppargs=["-O3"])
    except RuntimeError as r:
        m = re.search("In instant.recompile: The module did not compile with command 'make VERBOSE=1', see '(.+?)'", r.message)
        if m:
            print m.group(1)
            os.system('cat %s' %m.group(1))
            exit()
    
    return compiled_module
compiled_module = _compile_code()

class UnassembledMatrix(compiled_module.UnassembledMatrix):
    def __init__(self, form, *args):
        from dolfin.fem.assembling import _create_dolfin_form
        compiled_module.UnassembledMatrix.__init__(self, _create_dolfin_form(form), *args)

if __name__ == '__main__':
    from dolfin import *
    import sys
    from cbcpost.utils.utils import cbc_print
    
    from numpy import random, array, set_printoptions

    N = int(sys.argv[1])
    mesh = UnitCubeMesh(N,N,N)
    
    V = FunctionSpace(mesh, "CG", 1)
    u,v = TrialFunction(V), TestFunction(V)
    
    print "Number of cells: ", mesh.num_cells()
    print "Function space dimensions: ", V.dim()

    U = []
    for i in range(mesh.geometry().dim()):
        U.append(Function(V))
    
    for i, _U in enumerate(U):
        _U.vector()[:] = random.random(V.dim())
        #_U.vector()[:] = array(range(V.dim()))
        #_U.vector()[:] = 1.0
    #U = [U1, U2, U3]
    
    form = inner(v, dot(as_vector(U), nabla_grad(u)))*dx
    
    tic()
    UM = UnassembledMatrix(form)
    print "Time to init: ", toc()

    M = Matrix()
    M2 = Matrix()
    
    avg = 0
    avgT2 = 0
    N = 20
    
    for i in xrange(N):
        cbc_print("*"*40)
        tic()
        #M2 = assemble(form)
        assemble(form, tensor=M)
        T1 = toc()
        cbc_print("Time spent ordinary assemble: %f"  %T1)

        tic()
        #UM.assemble(M2, U)
        UM.assemble(M2)
        T2 = toc()
        cbc_print("Time spent new assemble: %f" %T2)
        cbc_print("Time ratio: %f" %(T2/T1))
        cbc_print("Errornorm: %g" %(M2-M).norm('linf'))
        #print (M2-M).norm('linf')==0.0
        
        cbc_print("Norm new: %f" %(M2.norm('frobenius')))
        #cbc_print("Memory usage: %f" %getMemoryUsage())
        avg += T2/T1
        avgT2 += T2
        #exit()
    
    avg = avg/N
    avgT2 = avgT2/N
    #T2 = toc()
    add_local_calls = V.dim()
    
    print "Average time ratio: ", avg
    print "Average time new assemble: ", avgT2

    from fenicstools import getMemoryUsage
    for i in range(1000):
        UM.assemble(M2, U)
        cbc_print("Memory usage: %f" %getMemoryUsage())
    