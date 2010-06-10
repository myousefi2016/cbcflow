__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2010-04-15"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# List of problems
problems = ["drivencavity", "channel", "taylorgreen", "cylinder", "beltrami", "aneurysm"]

# Wrapper problem classes
def Problem(name, options):
    "Return problem instance for given problem name"
    exec("from %s import Problem as NamedProblem" % name)
    return NamedProblem(options)
