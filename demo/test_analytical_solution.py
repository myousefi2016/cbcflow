
import sys; sys.path.insert(0,"../site-packages")

from headflow import *
import beltrami

params = beltrami.Beltrami.default_user_params()
analytical_solution = beltrami.Beltrami.analytical_solution

analyzer1 = AnalyticalSolutionAnalyzer()
analyzer1.params["saveparams"]["save"] = True
analyzer2 = EnergyAnalyzer()
analyzer2.params["saveparams"]["save"] = True


ppfield_pd = ParamDict(
    saveparams=ParamDict(
        save=True,
        ),
    timeparams=ParamDict(
        step_frequency=10,
        )
    )
analyzer3 = WSS(params=ppfield_pd)


print params

for N in [4, 8, 16]: 
  for dt in [0.1, 0.05, 0.025]: 
    params["N"] = N 
    params["dt"] = dt 

    pp = NSPostProcessor({"casedir":"results/beltrami_ipcs_segregated/N=%d/dt=%e" % (N,dt)})
    pp.add_field(analyzer1)
    pp.add_field(analyzer2)
    pp.add_field(analyzer3)

    p = beltrami.Beltrami(params)
    scheme = SegregatedIPCS(None)
    nssolver = NSSolver(p, scheme, pp)  
    nssolver.solve()

    pp = NSPostProcessor({"casedir":"results/beltrami_ipcs/N=%d/dt=%e" % (N,dt)})
    pp.add_field(analyzer1)
    pp.add_field(analyzer2)
    pp.add_field(analyzer3)

    p = beltrami.Beltrami(params)
    scheme = IPCS(None)
    nssolver = NSSolver(p, scheme, pp)  
    nssolver.solve()

     
  

