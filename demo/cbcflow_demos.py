#!/usr/bin/env python
import sys, os
import os.path
from os.path import split, join, abspath
from glob import glob
from cbcflow import NSProblem

parents = ("documented", "undocumented")

def find_demo_paths():
    paths = {}
    for p in parents:
        paths[p] = [path for path in glob(join(p, "*")) if os.path.isdir(path)]
    return paths

def find_problems(module):
    return [v for v in vars(module).values()
            if isinstance(v, type) and issubclass(v, NSProblem) and not v is NSProblem]

def print_demo_paths(paths):
    print "Documented demos:"
    print '\n'.join('    '+p for p in paths["documented"])
    print
    print "Undocumented demos:"
    print '\n'.join('    '+p for p in paths["undocumented"])
    print

def print_demo_info(fullname, module, problems):
    print "="*80
    print "In", fullname
    if module.__doc__ is not None:
        print module.__doc__
    for Problem in problems:
        print
        print "Class %s:" % (Problem.__name__,)
        print Problem.__doc__
        print

def list_demos():
    demo_paths = find_demo_paths()
    print_demo_paths(demo_paths)

def list_demo_docs():
    demo_paths = find_demo_paths()
    for parent in parents:
        for path in demo_paths[parent]:
            # Find demo python name
            fullname, = glob(join(path, "*.py"))

            # Import demo and find problem classes
            sys.path.insert(0, path)
            module_name = split(fullname)[-1].replace(".py", "")
            module = __import__(module_name)
            problems = find_problems(module)
            print_demo_info(fullname, module, problems)

def run_demos():
    root = abspath(os.curdir)
    demo_paths = find_demo_paths()
    for parent in parents:
        for path in demo_paths[parent]:
            # Find demo python name
            fullname, = glob(join(path, "*.py"))
            pyname = split(fullname)[-1]

            # Execute demo
            demoroot = abspath(join(root, path))
            os.chdir(demoroot)
            try:
                cmd = "./%s > demo.log" % pyname
                print "In '%s', running '%s':" % (demoroot, cmd)
                #status, output = get_status_output(cmd) # FIXME: Execute demos! Maybe move this code to the tests directory.
            except:
                print "Executing %s failed!" % fullname
            os.chdir(root)

def main(args):
    if '--list' in args:
        list_demos()
    elif '--docs' in args:
        list_demo_docs()
    elif '--run' in args:
        run_demos()
    else:
        print "Please provide either --list, --docs or --run."

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
