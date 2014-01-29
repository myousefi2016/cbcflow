#!/usr/bin/env python
import sys
import os.path
from os.path import join
from glob import glob
from cbcflow import NSProblem

def main(args):
    parents = ("documented", "undocumented")
    paths = {}
    for p in parents:
        paths[p] = [path for path in glob(join(p, "*")) if os.path.isdir(path)]

    print
    print "Documented demos:"
    print '\n'.join(paths["documented"])
    print
    print "Undocumented demos:"
    print '\n'.join(paths["undocumented"])
    print

    for p in parents:
        for path in paths[p]:
            sys.path.insert(0, path)
            pyname, = glob(join(path, "*.py"))
            module_name = os.path.split(pyname)[-1].replace(".py", "")

            module = __import__(module_name)
            problems = [v for v in vars(module).values() if isinstance(v, NSProblem)]

            if problems:
                Problem, = problems
                print "In", path
                print "Class", k
                print Problem.__doc__

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
