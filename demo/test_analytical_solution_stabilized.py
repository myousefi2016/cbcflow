
import sys; sys.path.insert(0,"../site-packages")

from headflow import *
import beltrami

Ns = [2, 4, 8]
dts = [0.1, 0.05, 0.025]
schemes = [IPCS(None), IPCS_Stable(None), IPCS_Stabilized(None), SegregatedIPCS(None)] 
schemes = [IPCS(None), IPCS_Stable(None), IPCS_Stabilized(None), SegregatedIPCS(None)] 

ppfield_pd = ParamDict(
    saveparams=ParamDict(
	save=True,
	),
    timeparams=ParamDict(
	step_frequency=10,
	)
    )
analyzer3 = WSS(params=ppfield_pd)

for mu in [1.0, 1.0e-3]: 
    params = beltrami.Beltrami.default_user_params()
    params["mu"] = mu 

    analytical_solution = beltrami.Beltrami.analytical_solution

    for scheme in schemes: 
	error_data = {}
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

