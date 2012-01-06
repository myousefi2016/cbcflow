#!/usr/bin/env python


# ... Template for command and casename
params = ['stationary', 'test_case', 'refinement_level',
          'boundary_layers', 'dt', 'max_t',]

template = './ns Challenge ipcs_opt \
segregated=True \
save_solution=True \
save_frequency=%(save_frequency)d \
save_xml=True \
store_probes=False \
plot_probes=False \
stationary=%(stationary)s \
test_case=%(test_case)d \
refinement_level=%(refinement_level)d \
boundary_layers=%(boundary_layers)s \
dt=%(dt)g \
max_t=%(max_t)g \
casename=%(casename)r'


fenics_stable_setup = '''
echo Using fenics setup from
ls -lst ${FENICS_STABLE_DIR}/share/dolfin/dolfin.conf
source ${FENICS_STABLE_DIR}/share/dolfin/dolfin.conf
fenics-versions
'''


# ... Utility code
jobs = []
names = []
def name(**values):
    return 'challenge__' + '__'.join('%s__%s' % (k,values[k]) for k in params)
def emit_job(**values):
    casename = name(**values)
    values['casename'] = casename
    jobs.append(template % values)
    names.append(casename)

debug = False

# ... Parameter sweep
smax_t = 0.3
pmax_t = 2.5
#pmax_t = 1.5
# Set dt to None to use computed cfl condition
#for dt in (1e-4,):
#for dt in (1e-3,1e-4,1e-5):
for dt in (2e-4,2e-5):
    save_frequency = 10 if dt is None else max(1,int(0.5+1.0/(400*dt)))

    #smax_t = pmax_t = dt*10 # Short timespan for debugging!

    #for refinement_level in (0,1,2):
    #for refinement_level in (0,):
    for refinement_level in (0,1,):
        #for boundary_layers in (False,True):
        for boundary_layers in (True,):

            stationary = True
            max_t = smax_t
            #for test_case in (1,):
            for test_case in (1,2,3,4):
                emit_job(**locals())
                if debug: break
            if debug: break

            stationary = False
            max_t = pmax_t
            #for test_case in (1,2,):
            for test_case in (1,2,):
                emit_job(**locals())
        if debug: break
    if debug: break

# ... Execute jobs!
if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print "Pass number of processors as first argument. Job list:"
        print '\n'.join(jobs)
    else:
        from dolfin_utils.pjobs import submit
        n = int(sys.argv[1]) 
        if n < 1:
            submit(jobs, name=names, serial=True)
        else:
            setup = ''#fenics_stable_setup # TODO: (msa) Test this now that correct .bash_profile has been configured
            prefix = "pmpirun.openmpi -n %d " % n
            jobs = [prefix + job for job in jobs]
            nn = (n+7)//8
            assert n <= nn * 8
            submit(jobs, name=names, nodes=(n+7)//8, ppn=8, walltime=24*7, setup=setup, keep_environment=False, dryrun=False)

