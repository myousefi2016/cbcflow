Things to do in headflow
************************

- The WIP (work in progress) lists is what we work on short term

- The backlog contains todos nobody are working on yet


WIP Martin
==========


WIP Oyvind
==========

- Add your own focus points here!


WIP Kent
========

- Add your own focus points here!


Backlog
=======

Known bugs
----------

- CoupledNonLinear has incorrect boundary terms (try channel demo to see)


General
-------

- Find a better name than headflow. Suggestions:

    - bioflow
    - fenflow
    - dolflow
    - dolflo

- Make an email list (using google groups?)

- Make a team and project repository on bitbucket

- Get buildbot running

- Add licence headers

- Rename schemes and their files to follow at least some resemblance
  of a convention.

- Reorganize submodule structure a bit?

- Document schemes with LaTeX and rst in docstrings

- Document demo problems with LaTeX and rst in docstrings

- Make a release


Parameters
----------

- Include and improve upon parameter sweep scripts. Cover this well in unit tests!


BCs
---

- Update boundary handling to fenics dev? Need to agree on a time to do this.

- Improve boundary condition utilities further (not sure what is needed yet).
  Cover these well in unit tests!


Postprocessing
--------------

- Add scripts for inspecting cases

- Handle transient dependencies

- Improve postprocessing further (not sure what is needed yet).
  Cover this well in unit tests!


Tools
-----

- Improve show_problem to be more useful.
  Hard to unit test with all the plotting though.

- Add a check_problem to do something similar to show_problem but without plotting.
  Cover with unit tests demonstrating various possible user errors

- Transform Kents demo/test_*.py into something reusable (analyze_problem_foo?)


Demos
-----

- Update the last demos for new problem interface

- Fix and validate all demo problems for all schemes


Tests
-----

- Add NSSolver unit test with mock scheme and postprocessor,
  testing that the solver steps correctly.

- Improve unit test suite further (not sure what is needed yet).
  Bonus points: Make it runnable on cluster.

- Setup a fast regression test suite (trivial meshes, few timesteps).
  Bonus points: Make it runnable on cluster.

- Setup a slow regression test suite (realistic meshes).
  Bonus points: Make it runnable on cluster.

- Setup a validation test suite with analytical problems.
  Bonus points: Make it runnable on cluster.

- Setup a benchmark suite.
  Bonus points: Make it runnable on cluster.

