
# Hack until we place these tests in a proper framework
import sys; sys.path.insert(0, "../demo")

from headflow import *

from beltrami import Beltrami as Problem

params = Problem.default_params()

#Scheme = SegregatedIPCS
Scheme = IPCS_Stable

ppfield_pd1 = ParamDict(
    saveparams=ParamDict(
        save=True,
        ),
    )

ppfield_pd = ParamDict(
    saveparams=ParamDict(
        save=True,
        ),
    timeparams=ParamDict(
        step_frequency=10,
        )
    )
analyzer3 = WSS(ppfield_pd)


Ns = [2, 4, 8, 16]
dts = [0.1, 0.05, 0.025, 0.0125]

error_data = {}
for params.N in Ns:
    for params.dt in dts:

        analyzer1 = AnalyticalSolutionAnalyzer(ppfield_pd1)
        analyzer2 = EnergyAnalyzer(ppfield_pd1)

        casedir = "results/beltrami_ipcs_segregated/N=%d/dt=%e" % (N,dt)
        pp = NSPostProcessor({"casedir":casedir})
        pp.add_field(analyzer1)
        pp.add_field(analyzer2)
        pp.add_field(analyzer3)

        p = Problem(params)
        scheme = Scheme()
        nssolver = NSSolver(p, scheme, pp)
        nssolver.solve()

        error_data[(N, dt)] = analyzer1.get_data()
#        print error_data[(N, dt)]["data"]["u1"]


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
