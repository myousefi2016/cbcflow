Things to do in headflow
************************

- The WIP (work in progress) lists is what we work on short term

- The backlog contains todos nobody are working on yet


WIP Martin
==========

- Add NSScheme unit test with mock problem stepping through each scheme for 2
  timesteps and checking types, using headflow.all_schemes and their default parameters.

- Consolidate NSProblem and NSProblem2 and update all schemes to follow:

  - Use Spaces in all schemes

  - Use assign*ics in all schemes

  - Call update at timestep 0 in all schemes

  - Use make*bcs in all schemes

  - Remove fetch_bcs from scheme class

  - Return local namespace from all schemes


WIP Oyvind
==========

- Add your own focus points here!


WIP Kent
========

- Add your own focus points here!


Backlog
=======

General
-------

- Rename schemes and their files to follow at least some resemblance
  of a convention.

- Add licence headers

- Get buildbot running

- Find a better name than headflow

- Make an email list (using google groups?)

- Make a team and project repository on bitbucket

- Make a release


Parameters
----------

- Add locking of parameter names to ParamDict, a serious correctness liability.
  Cover with unit tests, making sure we never use the wrong parameter names!

- Improve Parameterized system: doesn't work that well with
  subclassing, and maybe don't mix subclass params into base namespace?

- Include and improve upon parameter sweep scripts.
  Cover this well in unit tests!


BCs
---

- Update boundary handling to fenics dev? Need to agree on a time to do this.

- Improve boundary condition utilities further (not sure what is needed yet).
  Cover these well in unit tests!


Postprocessing
--------------

- Improve postprocessing further (not sure what is needed yet).
  Cover this well in unit tests!

- Move data flow code from ppfield baseclass to nspostprocessor

- Handle transient dependencies

- Add scripts for inspecting cases


Tools
-----

- Improve show_problem to be more useful.
  Hard to unit test with all the plotting though.

- Add a check_problem to do something similar to show_problem but without plotting.
  Cover with unit tests demonstrating various possible user errors



Tests
-----

- Add NSSolver unit test with mock scheme and postprocessor,
  testing that the solver steps correctly.

- Improve unit test suite further (not sure what is needed yet).
  Bonus points: Make it runnable on cluster.

- Setup a validation test suite with analytical problems.
  Bonus points: Make it runnable on cluster.

- Setup a fast regression test suite (trivial meshes, few timesteps).
  Bonus points: Make it runnable on cluster.

- Setup a slow regression test suite (realistic meshes).
  Bonus points: Make it runnable on cluster.

- Setup a benchmark suite.
  Bonus points: Make it runnable on cluster.

- Setup a benchmark suite.
  Bonus points: Make it runnable on cluster.

- Make all suites run on buildbot!
