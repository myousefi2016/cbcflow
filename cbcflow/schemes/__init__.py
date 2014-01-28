"""A collection of Navier-Stokes schemes."""

from .official import official_schemes
for f in official_schemes:
    exec("from .official import %s" % (f,))

from .experimental import experimental_schemes
for f in experimental_schemes:
    exec("from .experimental import %s" % (f,))

all_schemes = official_schemes + experimental_schemes

def show_schemes():
    "Lists which schemes are available."
    print "Official schemes available:"
    print "\n".join("    " + f for f in official_schemes)
    print "Experimental schemes available:"
    print "\n".join("    " + f for f in experimental_schemes)

__all__ = (
    ["show_schemes", "all_schemes", "official_schemes", "experimental_schemes"] +
    all_schemes
    )
