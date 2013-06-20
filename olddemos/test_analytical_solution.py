
from headflow import *
import beltrami

params = beltrami.Beltrami.default_user_params()
analytical_solution = beltrami.Beltrami.analytical_solution


ppfield_pd = ParamDict(
    saveparams=ParamDict(
        save=True,
        ),
    timeparams=ParamDict(
        step_frequency=10,
        )
    )
analyzer3 = WSS(params=ppfield_pd)


error_data = {}
Ns = [2, 4, 8, 16]
dts = [0.1, 0.05, 0.025, 0.0125]
for N in Ns:
    for dt in dts:
        params["N"] = N
        params["dt"] = dt

        analyzer1 = AnalyticalSolutionAnalyzer()
        analyzer1.params["saveparams"]["save"] = True
        analyzer2 = EnergyAnalyzer()
        analyzer2.params["saveparams"]["save"] = True

        pp = NSPostProcessor({"casedir":"results/beltrami_ipcs_segregated/N=%d/dt=%e" % (N,dt)})
        pp.add_field(analyzer1)
        pp.add_field(analyzer2)
        pp.add_field(analyzer3)

        p = beltrami.Beltrami(params)
#        scheme = SegregatedIPCS()
        scheme = IPCS_Stable()
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
