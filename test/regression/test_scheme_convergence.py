
import pytest

import logging
import os, sys, shelve
from hashlib import sha1
from numpy import sqrt

from ufl.classes import ListTensor
from dolfin import *

from cbcflow import *

import time

set_log_level(100)

class NoOutput():
    def write(self, s):
        pass
    
    def flush(self):
        pass

def l2norm(val, ref):
    if isinstance(val, (float, int)):
        err = (val-ref)**2
        ref = ref**2

    elif isinstance(val, Function):
        if isinstance(ref, list):
            ref = as_vector(ref)
        #degree = val.ufl_element().degree()
        # TODO: Ensure higher degree for ref here, or leave it up the the problem
        err = assemble((val-ref)**2*dx)
        ref = assemble(ref**2*dx, mesh=val.function_space().mesh())

    return sqrt(err/ref) if ref > 1e-14 else sqrt(err)

def write_output_data(filename, metadata, values):
    path = os.path.abspath(os.path.join(os.path.split(__file__)[0], "..", "output"))
    ref_filename = os.path.join(path, filename)
    if not os.path.isdir(path):
        os.mkdir(path)

    f = shelve.open(ref_filename, 'n')
    f["metadata"] = metadata
    f["values"] = values
    f.close()

def read_reference_data(filename):
    path = os.path.abspath(os.path.join(os.path.split(__file__)[0], "..", "cbcflow-reference-data"))
    ref_filename = os.path.join(path, filename)
    if os.path.isfile(ref_filename):
        ref = shelve.open(ref_filename, 'r')
        return ref["metadata"], ref["values"]
    else:
        return {}, {}

def test_run_problem(problem_factory, scheme_factory, refinement_level, dt):

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
    values = {f.name: 1e16 for f in test_fields}
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
        num_dofs = spaces.V.dim()+spaces.Q.dim()

        references = problem.test_references(spaces, t)
        if references:
            assert len(references) == len(test_fields)
            for field, ref in zip(test_fields, references):
                value = pp.get(field.name)
                values[field.name] = l2norm(value, ref)
        else:
            for field in test_fields:
                value = float(pp.get(field.name)) # Only support scalar values in reference data
                values[field.name] = value

    except RuntimeError as re:
        print re.message

    # Enable printing again, and print values
    sys.stdout = original_stdout

    print "** dofs: %d Time spent: %f" % (num_dofs , T)

    #assert values, "No values calculated. Solver most likely failed."
    if all([v==1e16 for v in values.values()]):
        print "No values calculated. Solver most likely failed."

    for tfname, err in values.items():
        print "**** Fieldname: %20s ** Error: %.8e" % (tfname, err)

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
    filename = hash.hexdigest() + ".db"

    # Always write to output
    write_output_data(filename, metadata, values)

    # Read reference data values
    ref_metadata, ref_values = read_reference_data(filename)
    assert ref_values != {}, "Found no reference data!"

    # Check each value against reference
    for key in values:
        if key in ref_values:

            # Compute absolute and relative errors
            abs_err = abs(values[key] - ref_values[key])
            if abs(ref_values[key]) > 1e-12:
                err = abs_err / abs(ref_values[key])
            else:
                err = abs_err

            # TODO: Find necessary condition of this check!

            # This one should be chosen such that it always passes when nothing has changed
            strict_tolerance = 1e-8
            if err > strict_tolerance:
                msg = "Error not matching reference with tolerance %e:\n    key=%s,  error=%e,  ref_error=%e  diff=%e,  relative=%e" % (
                    strict_tolerance, key, values[key], ref_values[key], abs_err, err)
                print msg

            # This one should be chosen so that it passes when we're happy
            loose_tolerance = 1e-3
            assert err < loose_tolerance

    # After comparing what we can, check that we have references for everything we computed
    assert set(values.keys()) == set(ref_values.keys()), "Value keys computed and in references are different."
