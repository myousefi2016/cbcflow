import shelve
import os
from hashlib import sha1
import sys
import pprint

def print_table(scheme_data, problem_data, data):

    # Print table header
    print "Scheme: ", scheme_data["name"]
    print "Parameters:"
    print scheme_data["params"]
    print
    print "Problem: ", problem_data["name"]
    print "Problem parameters:"
    print problem_data["params"]
    
    fields = data.pop("fields")

    # Print table (cw=column width)    
    cw = 12
    for f in fields:
        print
        print f
        for i, dt in enumerate(reversed(sorted(data.keys()))):
            if i==0:
                for refinement_level in [""]+sorted(data[dt].keys()):
                    sys.stdout.write('{0:{align}{w}}'.format(refinement_level, w=cw, align="^"))
                sys.stdout.write('\n')
            sys.stdout.write('{0:{align}{w}}'.format(dt, w=cw, align="^"))
            for refinement_level in sorted(data[dt].keys()):
                if data[dt][refinement_level] and f in data[dt][refinement_level]:
                    err = "%.4e" %data[dt][refinement_level][f]
                else:
                    err = "N/A"
                
                sys.stdout.write('{0:{align}{w}}'.format(err, w=cw, align="^"))
            sys.stdout.write('\n')
    print
    print "-"*80
    print
    
    


L = []

if __name__ == '__main__':
    
    table_dict = {}
    
    # Read all files from reference folder, and sort into table dict
    folder = 'reference'
    for f in os.listdir(folder):
        datafile = dict(shelve.open(os.path.join(folder, f), 'r'))
        
        problem_data = datafile["metadata"]["problem"]
        scheme_data = datafile["metadata"]["scheme"]
        
        refinement_level = problem_data["params"].pop("refinement_level")
        dt = problem_data["params"].pop("dt")
        
        scheme_hash = sha1(str(scheme_data)).hexdigest()
        problem_hash = sha1(str(datafile["metadata"]["problem"])).hexdigest()
        
        num_dofs = datafile["metadata"]["num_dofs"]
        time = datafile["metadata"]["time"]
        
        if not scheme_hash in table_dict:
            table_dict[scheme_hash] = {}
            table_dict[scheme_hash]["scheme"] = datafile["metadata"]["scheme"]
        
        if not problem_hash in table_dict[scheme_hash]:
            table_dict[scheme_hash][problem_hash] = {}
            table_dict[scheme_hash][problem_hash]["problem"] = problem_data
            table_dict[scheme_hash][problem_hash]["fields"] = []
        
        if dt not in table_dict[scheme_hash][problem_hash]:
            table_dict[scheme_hash][problem_hash][dt] = {}
        if refinement_level not in table_dict[scheme_hash][problem_hash][dt]:
            table_dict[scheme_hash][problem_hash][dt][refinement_level] = {}

        table_dict[scheme_hash][problem_hash][dt][refinement_level] = datafile["errors"]
        if datafile["errors"]:
            table_dict[scheme_hash][problem_hash]["fields"] = list(set(table_dict[scheme_hash][problem_hash]["fields"]+datafile["errors"].keys()))
    
    # Print all tables
    for scheme_hash in table_dict.keys():
        scheme_data = table_dict[scheme_hash].pop("scheme")
        for problem_hash in table_dict[scheme_hash]:
            problem_data = table_dict[scheme_hash][problem_hash].pop("problem")
            print_table(scheme_data, problem_data, table_dict[scheme_hash][problem_hash])
