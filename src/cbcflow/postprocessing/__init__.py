
# These are all the built-in PPField classes
all_fields = [
    # The basic solution fields:
    "Velocity", "Pressure", "PhysicalPressure",
    # The basic solution derivative fields:
    "VelocityGradient", "VelocityCurl", "VelocityDivergence",
    "FlowRate",
    "PressureGradient",
    # Material model specific fields:
    "Strain", "Stress",
    # Derived fields on boundary:
    "WSS",
    # TODO: More indicators: "OSI", ...
    # Fields used for inspecting parts of fields
    "SubFunction", "PointEval",
    # Fields used for ???:
    "LocalCfl",
    "Q", "Delta", "Lambda2",
    # Errors w.r.t. analytical solutions, used for testing:
    "AnalyticalVelocity", "AnalyticalPressure",
    "VelocityError", "PressureError",
    # Not sure exactly what this is for? Used in tests
    "KineticEnergy",
    "EnergyAnalyzer",
    # Metafields, applied to other fields:
    "TimeIntegral", "TimeDerivative", "SecondTimeDerivative",
    "L2norm", "H1norm", "H1seminorm", "Linfnorm",
    "DiffL2norm", "DiffH1norm", "DiffH1seminorm",
    "RunningAvg", "RunningMin", "RunningMax", "RunningL2norm",
    "DomainAvg", "BoundaryAvg", "Magnitude",
    # TODO: More metafields:
    #"Gradient", "Divergence", "Curl",
    ]

def show_fields():
    print "Postprocessing fields available:"
    print "\n".join("    " + f for f in all_fields)

# Import field classes from modules with same name
for f in all_fields:
    exec("from .%s import %s" %(f,f))
from .PPField import PPField
from .MetaPPField import MetaPPField, MetaPPField2

# Make a mapping from name to type, for use in NSPostProcessor
field_classes = { f: eval(f) for f in all_fields }
assert all(issubclass(c, PPField) for c in field_classes.itervalues())

# Define symbols to export
__all__ = ["PPField", "MetaPPField", "MetaPPField2", "field_classes", "all_fields"] + all_fields
