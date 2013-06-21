
# Hack until we place these tests in a proper framework
import sys; sys.path.insert(0, "../demo")

from headflow import *
from math import sqrt
import dolfin

from flow_around_cylinder import FlowAroundACylinder as Problem

dolfin.parameters["allow_extrapolation"] = True


#Ns = [2, 4]
Ns = [32, 64]
#dts = [0.1, 0.05]
dts = [0.025]

schemes = [IPCS(), IPCS_Stable(), IPCS_Stabilized(), SegregatedIPCS()]
schemes = [IPCS_Stabilized(), IPCS_Stable()]

ppfield_pd = ParamDict(
    saveparams=ParamDict(
        save=True,
        ),
    timeparams=ParamDict(
        step_frequency=10,
        )
    )

mus = [1.0e-2]
for mu in mus:
    params = Problem.default_params()
    params["mu"] = mu

    for scheme in schemes:
        velocity = {}
        for N in Ns:
            for dt in dts:
#                try:
                if 1:
                    params["N"] = N
                    params["dt"] = dt

                    p = Problem(params)

                    analyzer1 = Velocity(params=ppfield_pd)
                    analyzer2 = QDeltaLambda2(params=ppfield_pd)

                    pp = NSPostProcessor({"casedir":"results/%s/%s/mu=%s/N=%d/dt=%e" % (p.shortname(), scheme.shortname(), str(mu), N,dt)})
                    pp.add_field(analyzer1)
                    pp.add_field(analyzer2)

                    nssolver = NSSolver(p, scheme, pp)
                    nssolver.solve()

        print ""
        print "scheme ", scheme, " mu ", mu

        print velocity.keys()
        for i in range(len(Ns[1:])):
            for dt in dts:
                try:
                    u_fine = velocity[(Ns[i+1],dt)]["data"]
                    u_coarse = velocity[(Ns[i], dt)]["data"]
                    uc = dolfin.interpolate(u_coarse, u_fine.function_space())
                    difference = dolfin.assemble(dolfin.inner(u_fine - uc, u_fine - uc)*dolfin.dx())
                    u_norm = dolfin.assemble(dolfin.inner(u_fine, u_fine)*dolfin.dx())
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
                    difference = dolfin.assemble(dolfin.inner(u_fine - uc, u_fine - uc)*dolfin.dx())
                    u_norm = dolfin.assemble(dolfin.inner(u_fine, u_fine)*dolfin.dx())
                    print "difference between dt ", dts[i+1] , " and ", dts[i], " at level ", N, " is ", difference, " versus u_norm ", u_norm
                except Exception as e:
                    print "Not able to compare ", N, dt
                    print e
