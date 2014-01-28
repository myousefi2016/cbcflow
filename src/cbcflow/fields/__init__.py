
# Base classes for fields
from .bases.PPField import PPField
from .bases.MetaPPField import MetaPPField
from .bases.MetaPPField2 import MetaPPField2

# Lists of available field names
from .basic import basic_fields
from .meta import meta_fields
all_fields = basic_fields + meta_fields

# Import field classes from modules with same name
for f in basic_fields:
    exec("from .basic.%s import %s" % (f, f))
for f in meta_fields:
    exec("from .meta.%s import %s" % (f, f))

# Make a mapping from name to type, for use in NSPostProcessor
field_classes = { f: eval(f) for f in basic_fields }
assert all(issubclass(c, PPField) for c in field_classes.itervalues())

def show_fields():
    "Lists which fields are available."
    print "Postprocessing fields available by name:"
    print "\n".join("    " + f for f in basic_fields)
    print "Postprocessing fields available with parameters:"
    print "\n".join("    " + f for f in meta_fields)

# Define symbols to export
__all__ = (
      ["PPField", "MetaPPField", "MetaPPField2"]
    + ["show_fields", "field_classes", "basic_fields", "meta_fields", "all_fields"]
    + all_fields
    )

