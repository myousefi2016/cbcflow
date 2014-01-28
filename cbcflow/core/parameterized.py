from __future__ import division
__author__ = "Martin Alnaes <martinal@simula.no>"
__date__ = "2013-04-26"
__copyright__ = "Copyright (C) 2013-2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from .paramdict import ParamDict

class Parameterized(object):
    "Core functionality for parameterized subclassable components."
    def __init__(self, params):
        self.params = self.default_params()
        self.params.replace(params)

        # Assert for each subclass that we have all keys,
        # i.e. no default_params functions have been skipped
        # in the inheritance chain
        ps = set(self.params.keys())
        for cl in type(self).mro()[:-2]: # Skip object and Parameterized
            assert len(set(cl.default_params().keys()) - ps) == 0

    # --- Default parameter functions ---

    @classmethod
    def default_params(cls):
        "Merges base and user params into one ParamDict."
        raise NotImplementedError("Missing default_params implementation for class %s" % (cls,))

    # --- Name functions ---

    @classmethod
    def shortname(cls):
        """Get a one-word description of what the class represents.

        By default uses class name."""
        return cls.__name__

    @classmethod
    def description(cls):
        """Get a one-sentence description of what the class represents.

        By default uses first line of class docstring."""
        d = cls.__doc__
        if d is None:
            return "Missing description."
        else:
            return d.split('\n')[0]

    def __str__(self):
        return "%s: %s" % (self.shortname(), self.description())
