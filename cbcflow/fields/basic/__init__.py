
# Fields that can be constructed just by name
basic_fields = [
    # The basic solution fields:
    "Velocity",
    "Pressure",
    "PhysicalPressure",

    # Errors w.r.t. analytical solution if provided by problem:
    "AnalyticalVelocity",
    "AnalyticalPressure",
    "VelocityError",
    "PressureError",

    # Derived fields:
    "VelocityGradient",
    "VelocityCurl",
    "VelocityDivergence",
    "StreamFunction",
    "PressureGradient",
    "Strain",
    "Stress",
    "WSS",
    "LocalCfl",
    "KineticEnergy",
    "Q",
    "Delta",
    "Lambda2",
    ]

