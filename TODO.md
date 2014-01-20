Things to do in cbcflow
************************

- The WIP (work in progress) lists is what we work on short term

- The backlog contains todos nobody are working on yet


Backlog
=======

General
-------

- Make an email list (using google groups?)

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

- Make utilities to convert peak velocities to flow rates under assumption of Womersley flow (using Jarles formulas)


Postprocessing
--------------

- Extract a Plotter and Storage classes from NSPostProcessor to separate concerns better

- Test multiple levels of time dependencies (e.g. TimeDerivative(TimeDerivative(foo)) == SecondTimeDerivative(foo))

- Improve postprocessing further (not sure what is needed yet).
  Cover this well in unit tests!


Tools
-----

- Transform Kents demo/test_*.py into something reusable (analyze_problem_foo?)


Demos
-----

- Add some more interesting demo problems


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

