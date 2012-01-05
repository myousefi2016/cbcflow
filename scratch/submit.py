#!/usr/bin/env python

# TODO: Implement save_interval as alternative to save_frequency
# TODO: Implement restarting
# TODO: Implement probing utils

# Template for command and casename
params = ['stationary', 'test_case', 'refinement_level', 'boundary_layers', 'dt', 'max_t']
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


# Utility code
jobs = []
names = []
def name(**values):
    return 'challenge__' + '__'.join('%s__%s' % (k,values[k]) for k in params)
def emit_job(**values):
    casename = name(**values)
    values['casename'] = casename
    jobs.append(template % values)
    names.append(casename)


# Parameter sweep
max_t = 0.0001 # 3.0
for dt in (1e-5,):#(1e-3,1e-4,1e-5): # None to use cfl condition
    save_frequency = 10 if dt is None else max(1,int(0.5+1.0/(400*dt)))

    for refinement_level in (0,):#(0,1,2):
        for boundary_layers in (False,):#(False,True):

            stationary = True
            for test_case in (1,):#(1,2,3,4):
                emit_job(**locals())

            stationary = False
            for test_case in (1,):#(1,2,):
                emit_job(**locals())


# Execute jobs
if __name__ == '__main__':
    print '\n'.join(jobs)
    print '\n'.join(names)
    # FIXME:
    #from dolfin_utils.pjobs import submit
    #submit(jobs, names=names, walltime=24*7, nodes=2, ppn=8, keep_environment=True)

