#!/usr/bin/env python
"""
TODO: Make a demo of this with full usage example.
"""

from cbcflow import Field
import dolfin
import numpy

class CustomPostprocessingField(Field):
    def before_first_compute(self, pp, spaces, problem):
        # Initialize internal variables once in this optional function
        self._my_private_variable = 0
        self._my_dependency_name = "Velocity"

    def compute(self, pp, spaces, problem):
        # Get dependencies using one of these patterns:
        u = pp.get(self._my_dependency_name)
        p = pp.get("Pressure")

        # Do computations, update internal variables if needed:
        self._my_private_variable += 1

        # Return some value, the returned object will not be touched:
        return self._my_private_variable

    def after_last_compute(self, pp, spaces, problem):
        # Perform and return some final computation after all timesteps have completed:
        value = 1.0 / self._my_private_variable
        return value
