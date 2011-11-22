#!/usr/bin/env python

__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2010-05-20"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2010.
# Last changed: 2010-05-20

from pylab import *
import sys

# Results file
filename = "results/results.log"

# Tags for solvers
def tags(solver):
    t = {"Chorin": 'r-v',
         "CSS1"  : 'y-^',
         "CSS2"  : 'g-D',
         "GRPC"  : 'k-s',
         "G2"    : 'b-p',
         "IPCS"         : 'c-o',
         "IPCS_p1p1"    : 'r-^',
         "IPCS_p1p1_seg": 'y-D',
         "IPCS_opt"     : 'g-s',
         "IPCS_opt_seg" : 'k-p',
         }
    if solver in t:
        return t[solver]
    print 'Unknown solver "%s"'%solver
    return 'm-o'

# Extract data from file
data = {}
for lno, line in enumerate(line for line in open(filename, "r").read().split("\n") if "," in line):

    # Extract data
    date, problem, solver, num_dof, cputime, wct,func, dt_ref, error = line.split(",")
    if float(error) == 0:
        print "Ignoring "+filename+":"+str(lno+1)+" -- error not recorded"
        continue
    problem = problem.strip()
    solver = solver.strip()

    # If problem names are given on the command line, show only those
    if len(sys.argv) > 1 and problem not in sys.argv[1:]:
        continue

    # Save data
    if not problem in data:
        data[problem] = {}
    if not solver in data[problem]:
        data[problem][solver] = ([], [], [], [])
    data[problem][solver][0].append(int(num_dof))
    data[problem][solver][1].append(float(cputime))
    data[problem][solver][2].append(float(error))
    data[problem][solver][3].append(float(func))

# Plot data
plot_kwargs = dict(linewidth=2, ms=10, alpha=1.0)
for i, (problem, problem_data) in enumerate(data.items()):

    # Create new plot window
    figure(i)

    # Plot results for each solver
    for solver, solver_data in problem_data.items():

        # Get data
       	num_dofs, cputimes, errors, funcs = solver_data

        # Set tag
        tag = tags(solver)

        # Plot
        if problem[0:5] == 'Cylin':
            subplot(211)
            loglog(num_dofs, cputimes, tag, **plot_kwargs)
            ylabel("CPU time", fontsize=15, color='black')
            grid(True)

            subplot(212)
            subplot(223)
            loglog(num_dofs, errors, tag, **plot_kwargs)
            ylabel("Errors", fontsize=15, color='black')
            grid(True)
            xlabel("Degrees of freedom", fontsize=15, color='black')

            subplot(224)
            semilogx(num_dofs, funcs, tag, **plot_kwargs)
            ylabel("Functional", fontsize=15, color='black')
            grid(True)

        elif problem[0:5] == 'Aneur':
            subplot(211)
            loglog(num_dofs, cputimes, tag, **plot_kwargs)
            ylabel("CPU time", fontsize=15, color='black')
            grid(True)

            subplot(212)
            loglog(num_dofs, errors, tag, **plot_kwargs)
            ylabel("Functional", fontsize=15, color='black')
            grid(True)

        else:
            subplot(211)
            loglog(num_dofs, cputimes, tag, **plot_kwargs)
            ylabel("CPU time", fontsize=15, color='black')
            grid(True)

            subplot(212)
            loglog(num_dofs, errors, tag, **plot_kwargs)
            ylabel("Errors", fontsize=15, color='black')
            grid(True)

    xlabel("Degrees of freedom", fontsize=15, color='black')

    # Set title and legend
    subplot(211)
    title(problem, fontsize=15, color='black')
    legend(problem_data, loc=2)

# Extract scatter plot data for CPU time vs error
points = {}
for problem, problem_data in data.items():
    mean_cputime_smallest_mesh = 1
    mean_error_smallest_mesh   = 1
    for solver, solver_data in problem_data.items():
        num_dofs, cputimes, errors, funcs = solver_data
        mean_cputime_smallest_mesh *= cputimes[0]
        mean_error_smallest_mesh   *= errors[0]
    mean_cputime_smallest_mesh **= 1.0/len(problem_data)
    mean_error_smallest_mesh   **= 1.0/len(problem_data)

    # Get scaled values
    for solver, solver_data in problem_data.items():
        num_dofs, cputimes, errors, funcs = solver_data
        num_levels = len(num_dofs)
        if not solver in points:
            points[solver] = [[], []]
        points[solver][0] += [(cputimes[i]/mean_cputime_smallest_mesh) for i in range(num_levels)] + [None]
        points[solver][1] += [(errors[i]/mean_error_smallest_mesh)     for i in range(num_levels)] + [None]

# Create scatter plot of CPU time vs error
figure(len(data))
hold(True)
grid(True)
xlabel("Scaled CPU time")
ylabel("Scaled error")
for solver in points:
    tag = tags(solver)
    #c = tags(solver)[0]
    #scatter(points[solver][0], points[solver][1], c=c, s=150, alpha=0.75)
    loglog(points[solver][0], points[solver][1], tag, **plot_kwargs)
title("Solver performance (scaled by a per-problem constant)", fontsize=15, color='black')
legend([solver for solver in points], loc=1)

show()
