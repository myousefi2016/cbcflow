# Script for running benchmark on cluster.
# Important: Must run combination of solver/problem to compile
# before submitting jobs on a higher refinement level.

# Submission is paused for 5 seconds to prevent small jobs from accessing the result.log file simultaneously.

# Kristian Valen-Sendstad, 2009

import os
from dolfin_utils.pjobs import submit
from time import *
#from pjobs import submit
if ".." not in os.environ['PYTHONPATH']:
    os.environ['PYTHONPATH'] = "..:" + os.environ['PYTHONPATH']
jobs = []

# Solvers and problems
solvers = ["chorin", "css1", "css2", "ipcs", "g2", "grpc"]
problems = ["drivencavity"  , "channel", "periodic", "beltrami", "cylinder", "aneurysm"]

# Number of refinement levels
refinements = [0,1,2,3,4,5,6,7,8,9,10]

# See output on screen
to_screen = True

# Max number of runtime days
days = 10

for k in range(len(refinements)):
        for j in range(len(problems)):
            jobs = []
            for i in range(len(solvers)):
                if to_screen==True:
                    print "python ns " + problems[j], solvers[i], "refinement_level=" +str(refinements[k]), ">&", problems[j] + "_" + solvers[i] + "_" + str(refinements[k])
                    jobs.append(("python ns " + problems[j] + " "+ solvers[i] + " " + "refinement_level=" +str(refinements[k]) + " "+ "debug=True" + ">&" + " " + problems[j]  + solvers[i] + "_" + str(refinements[k])))

            print ''
            if problems[j] == "beltrami":
                submit(jobs, nodes=1, ppn=8,  keep_environment=True, walltime=24*days, dryrun=False)
            else:
                submit(jobs, nodes=1, ppn=2,  keep_environment=True, walltime=24*days, dryrun=False)
            print ''
            print 'Submitted jobs:' , len(jobs)
            print ''
            sleep(5)





# For larger jobs (More memory demanding jobs, might require all available memory on node =>
# args nodes=1, ppn=8 allocates entire node for one job.)
#submit(jobs, nodes=1 , ppn=8 , keep_environment=True, walltime=24*days, dryrun=False)
