#!/usr/bin/env python


# ... Template for command and casename
params = ['stationary', 'test_case', 'refinement_level',
          'boundary_layers', 'dt', 'max_t',]

template = './ns Challenge ipcs_opt \
segregated=True \
save_solution=True \
save_frequency=%(save_frequency)d \
store_probes=False \
plot_probes=False \
stationary=%(stationary)s \
test_case=%(test_case)d \
refinement_level=%(refinement_level)d \
boundary_layers=%(boundary_layers)s \
dt=%(dt)g \
max_t=%(max_t)g \
casename=%(casename)r'


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


# ... Parameter sweep
#max_t = 3.5
#max_t = 1.5
# Set dt to None to use computed cfl condition
#for dt in (1e-3,1e-4,1e-5):
for dt in (1e-4,):
    save_frequency = 10 if dt is None else max(1,int(0.5+1.0/(400*dt)))

    max_t = dt*10 # Short timespan for debugging!

    #for refinement_level in (0,1,2):
    for refinement_level in (0,):
        #for boundary_layers in (False,True):
        for boundary_layers in (False,):

            stationary = True
            #for test_case in (1,2,3,4):
            for test_case in (1,):
                emit_job(**locals())

            stationary = False
            #for test_case in (1,2,):
            for test_case in (1,):
                emit_job(**locals())


# ... Execute jobs!
if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print "Pass number of processors as first argument. Job list:"
        print '\n'.join(jobs)
    else:
        from dolfin_utils.pjobs import submit
        n = int(sys.argv[1]) 
        if n == 1:
            submit(jobs, names=names, serial=True)
        else:
            prefix = "pmpirun.openmpi -n %d " % n
            jobs = [prefix + job for job in jobs]
            nn = (n+7)//8
            assert n <= nn * 8
            submit(jobs, names=names, nodes=(n+7)//8, ppn=8, walltime=24*7, keep_environment=True)

