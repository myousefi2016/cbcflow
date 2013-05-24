
import sys; sys.path.insert(0,"../site-packages")

from headflow import *
import cylinder 
from math import sqrt
import dolfin

dolfin.parameters["allow_extrapolation"] = True


#Ns = [2, 4]
Ns = [2, 4, 8 , 16]
#dts = [0.1, 0.05]
dts = [0.1, 0.05, 0.025]

schemes = [IPCS(None), IPCS_Stable(None), IPCS_Stabilized(None), SegregatedIPCS(None)] 
schemes = [IPCS_Stabilized(None), IPCS(None)] 

ppfield_pd = ParamDict(
    saveparams=ParamDict(
	save=True,
	),
    timeparams=ParamDict(
	step_frequency=10,
	)
    )

mus = [1.0, 1.0e-3]
mus = [1.0e-3]
for mu in mus: 
    params = cylinder.FlowAroundACylinder.default_user_params()
    params["mu"] = mu 

    for scheme in schemes: 
	velocity = {}
	for N in Ns: 
	    for dt in dts: 
                try:  
		    params["N"] = N 
		    params["dt"] = dt 

		    p = cylinder.FlowAroundACylinder(params)

		    analyzer = Velocity(params=ppfield_pd)

		    pp = NSPostProcessor({"casedir":"results/%s/%s/N=%d/dt=%e" % (str(p), str(scheme), N,dt)})
		    pp.add_field(analyzer)

		    nssolver = NSSolver(p, scheme, pp)  
		    nssolver.solve()

		    velocity [(N, dt)] = analyzer.get_data()
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
		    difference = dolfin.assemble(dolfin.inner(u_fine - uc, u_fine - uc)*dolfin.dx())
		    print "difference between level ", N, " and ", N-1, " with dt ", dt, " is ", difference
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
		    print "difference between level ", N, " and ", N-1, " with dt ", dt, " is ", difference
        
		except Exception as e: 
                    print "Not able to compare ", N, dt
                    print e


