
# Hack until we place these tests in a proper framework
import sys; sys.path.insert(0, "../demo")

from headflow import *

from beltrami import Beltrami as Problem

#Ns = [2, 4]
Ns = [2, 4, 8]
#dts = [0.1, 0.05]
dts = [0.5]
dts = [0.1, 0.05]

schemes = [IPCS(), IPCS_Stable(), IPCS_Stabilized(), SegregatedIPCS()]
schemes = [IPCS_Stabilized(), IPCS() ]

ppfield_pd = ParamDict(
    saveparams=ParamDict(
        save=True,
        ),
    timeparams=ParamDict(
        step_frequency=10,
        )
    )
analyzer3 = WSS(params=ppfield_pd)

mus = [1.0, 1.0e-3]
mus = [1.0e-3]
for mu in mus:
    params = Problem.default_params()
    params["mu"] = mu

    for scheme in schemes:
        error_data = {}
        for N in Ns:
            for dt in dts:
                params["N"] = N
                params["dt"] = dt

                p = Problem(params)

                analyzer1 = AnalyticalSolutionAnalyzer()
                analyzer1.params["saveparams"]["save"] = True
                analyzer2 = EnergyAnalyzer()
                analyzer2.params["saveparams"]["save"] = True

                pp = NSPostProcessor({"casedir":"results/%s/%s/N=%d/dt=%e" % (p.shortname(), scheme.shortname(), N,dt)})
                pp.add_field(analyzer1)
                pp.add_field(analyzer2)
                pp.add_field(analyzer3)

                nssolver = NSSolver(p, scheme, pp)
                nssolver.solve()

                error_data[(N, dt)] = analyzer1.get_data()

        print ""
        print ""
        print "scheme ", scheme, " mu ", mu
        for x in ["u0", "u1", "u2", "p"]:
            print ""
            print ""
            print x
            print "dt ",
            for dt in dts:
                print dt,
            for N in Ns:
                print "\nN ", N,
                for dt in dts:
                    print " %2.2e " % error_data[(N, dt)]["data"][x],

        print ""
        print ""
