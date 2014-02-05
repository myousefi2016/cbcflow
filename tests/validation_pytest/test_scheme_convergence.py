from cbcflow import *
import logging
from numpy import sqrt
import sys
import os
import shelve
from hashlib import sha1
import time
from ufl.tensors import ListTensor
from dolfin import *
import pytest

set_log_level(100)

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
    def test_run_problem(self, problem_factory, scheme_factory, refinement_level, dt):

        problem = problem_factory(refinement_level, dt)
        scheme = scheme_factory()

        print 
        print "."*100
        print "** Problem: %16s ** Solver: %16s" % (problem.shortname(), scheme.shortname())
        print "** Refinement level: %2d ** dt: %.8f" % (refinement_level, dt)
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
        #sys.stdout = NoOutput()
        
        try:
            t1 = time.time()
            ns = solver.solve()
            t2 = time.time()
            T = t2-t1

            spaces = ns["spaces"]
            t = float(ns["t"])

            errors = {}
            references = problem.test_references(spaces, t)
            for field, ref in zip(test_fields, references):
                value = pp.get(field.name)

                # TODO: Use field
                errors[field.name] = l2norm(value, ref, spaces)

            num_dofs = spaces.V.dim()+spaces.Q.dim()

        except RuntimeError as re:
            print re.message
            pass

        # Enable printing again, and print errors
        sys.stdout = original_stdout

        print "** dofs: %d Time spent: %f" % (num_dofs , T)
        if errors:
            for tfname, err in errors.items():
                print "**** Fieldname: %20s ** Error: %.8e" %(tfname, err)
        else:
            print "**** No errors calculated. Solver most likely failed


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

        filename = hash.hexdigest()+".db"
        ref_filename = os.path.join("reference", filename)
        out_filename = os.path.join("output", filename)

        type = pytest.config.option.type
        if type == "update":
            # Write to reference/
            if not os.path.isdir("reference"):
                os.mkdir("reference")

            f = shelve.open(ref_filename, 'n')
            f["metadata"] = metadata
            f["errors"] = errors
            f.close()

        elif type == "changed":
            # Write to output/
            if not os.path.isdir("output"):
                os.mkdir("output")

            f = shelve.open(out_filename)
            f["metadata"] = metadata
            f["errors"] = errors
            f.close()

            # Check against reference
            assert os.path.isfile(ref_filename), "Unable to find file for case (hash=%s)" % hash.hexdigest()
            ref = shelve.open(ref_filename, 'r')
            ref_errors = ref["errors"]

            assert set(errors.keys()) == set(ref_errors.keys()), "Different set of errors computed."

            for key in errors:
                print abs(errors[key]-ref_errors[key])
                print abs(errors[key]-ref_errors[key])/abs(ref_errors[key])
                # TODO: Find necessary condition of this check!
                assert abs(errors[key]-ref_errors[key]) < 1e-6, "Error not matching reference: key=%s error=%e, ref_error=%e (diff=%e)" %(key, errors[key], ref_errors[key], abs(errors[key]-ref_errors[key]))
