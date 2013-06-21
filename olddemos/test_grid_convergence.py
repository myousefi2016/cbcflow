
# Hack until we place these tests in a proper framework
import sys; sys.path.insert(0, "../demo")

from headflow import *
from math import sqrt
import dolfin

from flow_around_cylinder import FlowAroundACylinder as Problem

dolfin.parameters["allow_extrapolation"] = True


mus = [0.5e-2]
#Ns = [2, 4]
Ns = [64, 128, 256]
#dts = [0.1, 0.05]
dts = [0.025, 0.0125]

schemes = [IPCS(), IPCS_Stable(), IPCS_Stabilized(), SegregatedIPCS()]
schemes = [IPCS_Stabilized({"theta":1.0}), IPCS_Stabilized({"theta":0.5}), IPCS_Stable()]

ppfield_pd = ParamDict(
    saveparams=ParamDict(
        save=True,
        ),
    timeparams=ParamDict(
        step_frequency=10,
        )
    )

params = Problem.default_params()
for params.mu in mus:

    for scheme in schemes:
        scheme_str = scheme.shortname()
        if isinstance(scheme, IPCS_Stabilized):
            scheme_str += str(scheme.params.theta)

        velocity = {}

        for params.N in Ns:
            for params.dt in dts:
                p = Problem(params)

                casedir = "results/%s/%s/mu=%s/N=%d/dt=%e" % (p.shortname(), scheme_str, str(mu), N,dt)
                pp = NSPostProcessor({"casedir":casedir})

                analyzer = Velocity(params=ppfield_pd)
                pp.add_field(analyzer)

                nssolver = NSSolver(p, scheme, pp)

                try:
                    nssolver.solve()
                    velocity [(N, dt, mu)] = analyzer.get_data()
                except Exception as e:
                    print "The scheme did not work "
                    print e

        print ""
        print "scheme ", scheme, " mu ", mu

        print velocity.keys()
        for i in range(len(Ns[1:])):
            for dt in dts:
                try:
                    u_fine = velocity[(Ns[i+1],dt)]["data"]
                    u_coarse = velocity[(Ns[i], dt)]["data"]
                    uc = dolfin.interpolate(u_coarse, u_fine.function_space())

                    difference = dolfin.assemble((u_fine - uc)**2*dolfin.dx())
                    u_norm = dolfin.assemble(u_fine**2*dolfin.dx())

                    print "difference between level ", Ns[i+1], " and ", Ns[i], " with dt ", dt, " is ", difference, " versus u_norm ", u_norm
                except Exception as e:
                    print "Not able to compare", N, dt
                    print e

        for i in range(len(dts[1:])):
            for N in Ns:
                try:
                    u_fine = velocity[(N,dts[i+1])]["data"]
                    u_coarse = velocity[(N, dts[i])]["data"]
                    uc = dolfin.interpolate(u_coarse, u_fine.function_space())

                    difference = dolfin.assemble((u_fine - uc)**2*dolfin.dx())
                    u_norm = dolfin.assemble(u_fine**2*dolfin.dx())

                    print "difference between dt ", dts[i+1] , " and ", dts[i], " at level ", N, " is ", difference, " versus u_norm ", u_norm

                except Exception as e:
                    print "Not able to compare ", N, dt
                    print e
