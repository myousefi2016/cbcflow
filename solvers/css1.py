__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2009-10-05"
__copyright__ = "Copyright (C) 2009-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from css import Solver as CSS

class Solver(CSS):
    "First-order consistent splitting scheme by Guermond and Shen."

    def __init__(self, options):
        CSS.__init__(self, options, 1)

    def __str__(self):
        return "CSS1"
