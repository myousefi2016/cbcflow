#!/usr/bin/env python

__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2010-05-20"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2010.
# Last changed: 2010-05-20

from pylab import *

# Results file
filename = "results/results.log"

# Tags for solvers
tags = {"Chorin": 'r-v', "CSS1": 'y-^', "CSS2": 'g-D', "IPCS": 'c-o', "GRPC": 'k-s', "G2": 'b-p'}

# Extract data from file
data = {}
for line in [line for line in open(filename, "r").read().split("\n") if "," in line]:

    # Extract data
    date, problem, solver, num_dof, cputime, wct,func, dt_ref, error = line.split(",")
    problem = problem.strip()
    solver = solver.strip()

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
for (i, problem) in []: #enumerate(data):

    # Create new plot window
    figure(i)

    # Plot results for each solver
    for solver in data[problem]:

        # Get data
       	num_dofs, cputimes, errors, funcs = data[problem][solver]

        # Set tag
        tag = tags[solver]

        # Plot
        if problem[0:5] == 'Cylin':
            subplot(211)
            loglog(num_dofs, cputimes, tag, linewidth=3, ms=10, alpha=1.0)
            ylabel("CPU time", fontsize=15, color='black')
            grid(True)

            subplot(212)
            subplot(223)
            plot(num_dofs, funcs, tag, linewidth=3, ms=10, alpha=1.0)
            ylabel("Functional", fontsize=15, color='black')
            grid(True)

            subplot(224)
            plot(num_dofs, funcs, tag, linewidth=3, ms=10, alpha=1.0)
            ylabel("Functional", fontsize=15, color='black')
            grid(True)

        elif problem[0:5] == 'Aneur':
            subplot(211)
            loglog(num_dofs, cputimes, tag, linewidth=3, ms=10, alpha=1.0)
            ylabel("CPU time", fontsize=15, color='black')
            grid(True)

            subplot(212)
            loglog(num_dofs, errors, tag, linewidth=3, ms=10, alpha=1.0)
            ylabel("Functional", fontsize=15, color='black')
            grid(True)

        else:
            subplot(211)
            loglog(num_dofs, cputimes, tag, linewidth=3, ms=10, alpha=1.0)
            ylabel("CPU time", fontsize=15, color='black')
            grid(True)

            subplot(212)
            loglog(num_dofs, errors, tag, linewidth=3, ms=10, alpha=1.0)
            ylabel("Errors", fontsize=15, color='black')
            grid(True)

            xlabel("Degrees of freedom", fontsize=15, color='black')

    # Set title and legend
    subplot(211)
    title(problem, fontsize=15, color='black')
    legend(data[problem], loc=2)

# Extract scatter plot data for CPU time vs error
points = {}
#data = {"Channel": data["Channel"]}
for problem in data:

    # Get maximum for scaling
    max_num_levels = max([len(data[problem][s][0]) for s in data[problem]])
    error_scale = [0.0 for i in range(max_num_levels)]
    cputime_scale = [0.0 for i in range(max_num_levels)]
    num_entries = [0 for i in range(max_num_levels)]
    for solver in data[problem]:
        num_dofs, cputimes, errors, funcs = data[problem][solver]
        for i in range(len(num_dofs)):
            error_scale[i] += errors[i] / num_dofs[i]
            cputime_scale[i] += cputimes[i] / num_dofs[i]
            num_entries[i] += 1
    error_scale = [error_scale[i] / num_entries[i] for i in range(max_num_levels)]
    cputime_scale = [cputime_scale[i] / num_entries[i] for i in range(max_num_levels)]

    # Get scaled values
    for solver in data[problem]:
        num_dofs, cputimes, errors, funcs = data[problem][solver]
        num_levels = len(num_dofs)
        if not solver in points:
            points[solver] = [[], []]
        points[solver][0] += [log(cputimes[i] / num_dofs[i] / cputime_scale[i]) for i in range(num_levels)]
        points[solver][1] += [log(errors[i] / num_dofs[i] / error_scale[i]) for i in range(num_levels)]

# Create scatter plot of CPU time vs error
figure(6)
hold(True)
grid(True)
xlabel("Scaled CPU time")
ylabel("Scaled error")
for solver in points:
    c = tags[solver][0]
    scatter(points[solver][0], points[solver][1], c=c, s=500.0, alpha=0.75)
title("Solver performance", fontsize=15, color='black')
legend([solver for solver in points], loc=3)

show()
