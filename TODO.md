

- Update boundary handling to fenics dev? Need to agree on a time to do this.


- Add NSSolver unit test with mock scheme and postprocessor,
  testing that the solver steps correctly.


- Rename schemes and their files to follow at least some resemblance
  of a convention.

- Consolidate NSProblem and NSProblem2 and update all schemes to follow.

- Add NSScheme unit test with mock problem stepping through each scheme for 2
  timesteps and checking types, using headflow.all_schemes and their default parameters.

- Use ics utilities in all schemes.

- Use make_bcs utilities in all schemes. Remove fetch_bcs from scheme class.

- Call update with timestep 0 before timeloop in all schemes.


- Add locking of parameter names to ParamDict, a serious correctness liability.
  Cover with unit tests, making sure we never use the wrong parameter names!

- Include and improve upon parameter sweep scripts.
  Cover this well in unit tests!

- Improve Parameterized system: doesn't work that well with
  subclassing, and maybe don't mix subclass params into base namespace?


- Improve show_problem to be more useful.
  Hard to unit test with all the plotting though.

- Add a check_problem to do something similar to show_problem but without plotting.
  Cover with unit tests demonstrating various possible user errors


- Improve postprocessing further (not sure what is needed yet).
  Cover this well in unit tests!


- Improve boundary condition utilities further (not sure what is needed yet).
  Cover these well in unit tests!


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


- Make suites run on buildbot!

- Make a release?

