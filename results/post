#!/usr/bin/env python

__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2010-05-20"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2010.
# Last changed: 2010-05-20

from pylab import *

# Results file
filename = "results.log"

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
for (i, problem) in enumerate(data):

    # Create new plot window
    figure(i)

    # Plot results for each solver
    for solver in data[problem]:

        # Get data
       	num_dofs, cputimes, errors, funcs = data[problem][solver]

        # Set tag
        tags = {"Chorin": 'r-v', "CSS1": 'y-^', "CSS2": 'g-D', "IPCS": 'c-o', "GRPC": 'k-s', "G2": 'b-p'}
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

show()
