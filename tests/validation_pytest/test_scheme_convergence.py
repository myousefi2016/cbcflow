from cbcflow import *
import logging
from numpy import sqrt
import sys
from dolfin import *
import os
import shelve
from hashlib import sha1
import time
from ufl.tensors import ListTensor

import pytest

class NoOutput():
    def write(self, s):
        pass
        
    
def l2norm(val, ref, spaces):
    if isinstance(val, (float, int)):
        err = (val-ref)**2
        ref = abs(ref)**2
    elif isinstance(val, Function):
        degree = val.ufl_element().degree()
        
        # Set degree as 3 higher than that if val if not already set
        if isinstance(ref, Expression):
            if not hasattr(ref, "degree"):
                ref.degree = degree+3
        elif isinstance(ref, list) or isinstance(ref, ListTensor):
            for ref_i in ref:
                if not hasattr(ref, "degree"):
                    ref_i.degree = degree+3
            ref = as_vector(ref)
        
        err = assemble(inner(val-ref,val-ref)*dx())
        ref = assemble(inner(ref, ref)*dx(), mesh=val.function_space().mesh())
    
    if ref > 1e-14:
        return sqrt(err/ref)
    else:
        return sqrt(err)
    
    
    
    

class TestConvergence():
    def test_run_problem(self,problem_factory, scheme_factory, refinement_level, dt):
        
        problem = problem_factory(refinement_level,dt)
        scheme = scheme_factory()
        
        print 
        print "."*100
        print "** Problem: %16s ** Solver: %16s" %(problem.shortname(), scheme.shortname())
        print "** Refinement level: %2d ** dt: %.8f" %(refinement_level, dt)
        print "**"
        
        pp = NSPostProcessor({"casedir": "test"})
        test_fields = problem.test_fields()
        pp.add_fields(test_fields)
        
        solver = NSSolver(problem, scheme, pp)
        
        # Define variables
        errors = None
        num_dofs = 0
        T = 0
        
        # Disable printing from solve
        # TODO: Move to NSSolver/set option
        original_stdout = sys.stdout
        sys.stdout = NoOutput()

        try:
            t1 = time.time()
            ns = solver.solve()
            t2 = time.time()
            T = t2-t1
            spaces = ns["spaces"]
            t = float(ns["t"])
            
            num_dofs = ns["spaces"].V.dim()+ns["spaces"].Q.dim()
            errors = {}
            references = problem.test_references(spaces, t)
            for i, tf in enumerate(test_fields):
                val = pp.get(tf.name)
                ref = references[i]
                
                # TODO: Use field
                errors[tf.name] = l2norm(val, ref, spaces)

        except RuntimeError as re:
            pass
        
        # Enable printing again, and print errors
        sys.stdout = original_stdout
        print "** dofs: %d Time spent: %f" %(num_dofs , T)
        if errors:
            for tfname, err in errors.items():
                print "**** Fieldname: %20s ** Error: %.8e" %(tfname, err)
        else:
            print "**** No errors calculated. Solver most likely failed."
        

        # Store solve metadata
        metadata = {}
        metadata["scheme"] = {}
        metadata["scheme"]["name"] = scheme.shortname()
        metadata["scheme"]["params"] = scheme.params
        metadata["problem"] = {}
        metadata["problem"]["name"] = problem.shortname()
        metadata["problem"]["params"] = problem.params
        metadata["num_dofs"] = num_dofs
        metadata["time"] = T
        
        # Find hash from problem and scheme name+parameters
        hash = sha1()
        hash.update(str(metadata["scheme"]))
        hash.update(str(metadata["problem"]))
        filename = os.path.join("reference", hash.hexdigest()+".db")
        
        type = pytest.config.option.type
        if type == "update":
            if not os.path.isdir("reference"):
                os.mkdir("reference")

            # Write to file
            f = shelve.open(filename, 'n')
            f["metadata"] = metadata
            f["errors"] = errors
            f.close()

        elif type == "changed":
            # Write to output
            if not os.path.isdir("output"):
                os.mkdir("output")
            f = shelve.open(os.path.join("output", hash.hexdigest()+".db"))
            f["metadata"] = metadata
            f["errors"] = errors
            f.close()
            
            # Check against reference
            assert os.path.isfile(filename), "Unable to find file for case (hash=%s)" %hash
            ref = shelve.open(filename, 'r')
            ref_errors = ref["errors"]
            if not errors:
                assert not ref_errors, "Unable to compute errors, but errors have been computed in reference"
            if errors:
                for key in errors:
                    assert key in ref_errors, "Unable to find key %s in reference" %key
                    assert abs(errors[key]-ref_errors[key]) < 1e-14, "Error not matching reference: key=%s error=%f, ref_error=%f" %(key, errors[key], ref_errors[key])
