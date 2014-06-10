#!/usr/bin/env python
import os
import sys
import glob
import shelve
from hashlib import sha1
from itertools import chain

def print_table(scheme_data, problem_data, case_data, write=sys.stdout.write):

    # Table entry formatting (cw=column width)
    cw = 12
    align = '{{0:{align}{w}}}'.format(w=cw, align='^')

    # Print table header
    print "Scheme: ", scheme_data["name"]
    print "Parameters:"
    print scheme_data["params"]
    print
    print "Problem: ", problem_data["name"]
    print "Problem parameters:"
    print problem_data["params"]
    
    # Loop over fields
    fields = case_data.pop("fields")
    for f in fields:
        print
        print f
        
        

        # Get row and column values
        dts = list(reversed(sorted(case_data.keys())))
        refinement_levels = sorted(set(chain(*(case_data[dt].keys() for dt in dts))))
        
        #print dts, refinement_levels

        # Write column headers
        if dts:
            for refinement_level in [""]+refinement_levels:
                write(align.format(refinement_level))
            write('\n')

        # Write the table rows
        for i, dt in enumerate(dts):
            # Row header
            write(align.format(dt))
            # Row values
            for refinement_level in refinement_levels:
                value = case_data[dt][refinement_level].get(f)
                err = "N/A" if value is None else ("%.4e"  % value)
                write(align.format(err))
            write('\n')
    print
    print "-"*80
    print

def read_tables(folder):
    # Read all files from reference folder, and sort into table dict
    table_dict = {}
    for f in glob.glob(os.path.join(folder, "*.db")):
        datafile = dict(shelve.open(f, 'r'))
        #import ipdb; ipdb.set_trace()

        # Unused:
        #num_dofs = datafile["metadata"]["num_dofs"]
        #time = datafile["metadata"]["time"]

        # Outer layer of table_dict is indexed by the hash of scheme data
        scheme_data = datafile["metadata"]["scheme"]
        scheme_hash = sha1(str(scheme_data)).hexdigest()
        if not scheme_hash in table_dict:
            table_dict[scheme_hash] = {}
            table_dict[scheme_hash]["scheme"] = scheme_data

        # The next layer of table_dict is indexed by the hash of problem data
        problem_data = datafile["metadata"]["problem"]
        dt = problem_data["params"].pop("dt")
        refinement_level = problem_data["params"].pop("refinement_level")
        problem_hash = sha1(str(problem_data)).hexdigest()
        if not problem_hash in table_dict[scheme_hash]:
            table_dict[scheme_hash][problem_hash] = {}
            table_dict[scheme_hash][problem_hash]["problem"] = problem_data
            table_dict[scheme_hash][problem_hash]["fields"] = []

        # This gives us a dict that will contain data for this scheme/problem
        case_data = table_dict[scheme_hash][problem_hash]

        # This will be indexed as case_data[dt][refinement_level] so build that structure
        if dt not in case_data:
            case_data[dt] = {}
        
        if refinement_level not in case_data[dt]:
            case_data[dt][refinement_level] = {}

        # Now store the test error values in this table entry
        case_data[dt][refinement_level] = datafile["values"]
        if datafile["values"]:
            case_data["fields"] = list(set(case_data["fields"] + datafile["values"].keys()))

    return table_dict

def print_all_tables(table_dict):
    "Print tables from table_dict[scheme_hash][problem_hash]. NB! Deletes values from table_dict."

    # For all schemes
    for scheme_hash in table_dict.keys():
        scheme_data = table_dict[scheme_hash].pop("scheme")

        # For all problems
        for problem_hash in table_dict[scheme_hash]:
            case_data = table_dict[scheme_hash][problem_hash]
            problem_data = case_data.pop("problem")

            # Print a table
            print_table(scheme_data, problem_data, case_data)

def main():
    #folder = 'cbcflow-reference-data'
    folder = 'output'
    table_dict = read_tables(folder)
    print_all_tables(table_dict)

if __name__ == '__main__':
    main()
