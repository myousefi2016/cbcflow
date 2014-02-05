
from cbcflow import *
from dolfin import *

from numpy import sqrt

import logging
import os, sys, shelve, time
from hashlib import sha1

from ufl.classes import ListTensor

import pytest

set_log_level(100)

class NoOutput():
    def write(self, s):
        pass

def l2norm(val, ref):
    if isinstance(val, (float, int)):
        err = (val-ref)**2
        ref = ref**2

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

        err = assemble((val-ref)**2*dx())
        ref = assemble(ref**2*dx(), mesh=val.function_space().mesh())

    return sqrt(err/ref) if ref > 1e-14 else sqrt(err)

def write_data(path, filename, metadata, errors):
    ref_filename = os.path.join(path, filename)
    if not os.path.isdir(path):
        os.mkdir(path)

    f = shelve.open(ref_filename, 'n')
    f["metadata"] = metadata
    f["errors"] = errors
    f.close()

def read_data(path, filename):
    ref_filename = os.path.join(path, filename)
    assert os.path.isfile(ref_filename), "Unable to find file %s in %s" % (filename, path)
    ref = shelve.open(ref_filename, 'r')
    return ref["metadata"], ref["errors"]

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
        errors = {}
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
            num_dofs = spaces.V.dim()+spaces.Q.dim()

            references = problem.test_references(spaces, t)
            for field, ref in zip(test_fields, references):
                value = pp.get(field.name)
                # TODO: Use field
                errors[field.name] = l2norm(value, ref)

        except RuntimeError as re:
            print re.message

        # Enable printing again, and print errors
        sys.stdout = original_stdout

        print "** dofs: %d Time spent: %f" % (num_dofs , T)
        for tfname, err in errors.items():
            print "**** Fieldname: %20s ** Error: %.8e" % (tfname, err)
        if not errors:
            print "**** No errors calculated. Solver most likely failed"

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
        # TODO: Why split errors out from metadata? More generic to just place it in there?
        #       We may want to reuse regression test code in more than one test function so generic would be nice.
        #metadata["errors"] = errors

        # Find hash from problem and scheme name+parameters
        hash = sha1()
        hash.update(str(metadata["scheme"]))
        hash.update(str(metadata["problem"]))

        filename = hash.hexdigest() + ".db"

        # FIXME: Drop this update action inside the test, it belongs in a script on the outside.
        # A typical workflow is this:
        # - Run tests. This step always writes to output, reads from reference, and reports the difference.
        # - If user wants to update references, call script to do so, which just copyies output/ over reference/, no rerun of tests necessary.

        type = pytest.config.option.type
        if type == "update":
            write_data("reference", filename, metadata, errors)

        elif type == "changed":
            write_data("output", filename, metadata, errors)
            ref_metadata, ref_errors = read_data("reference", filename)

            # Check against reference
            assert set(errors.keys()) == set(ref_errors.keys()), "Different set of errors computed."
            for key in errors:
                print abs(errors[key]-ref_errors[key])
                print abs(errors[key]-ref_errors[key])/abs(ref_errors[key])
                # TODO: Find necessary condition of this check!
                #       Suggestion: FFC has two tolerances, one strict and one acceptable. The acceptable one is enforced while the strict one is just reported. A bit like errors vs warnings.
                assert abs(errors[key]-ref_errors[key]) < 1e-6, "Error not matching reference: key=%s error=%e, ref_error=%e (diff=%e)" % (key, errors[key], ref_errors[key], abs(errors[key]-ref_errors[key]))
