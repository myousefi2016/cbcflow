
# Hack until we place these tests in a proper framework
import sys; sys.path.insert(0, "../demo") # FIXME: Assuming run from the headflow/tests/ directory!

from headflow import *
from math import sqrt
import dolfin
import unittest

dolfin.parameters["allow_extrapolation"] = True

def run_convergence_sweep(Scheme, scheme_params, scheme_str,
                          Problem, problem_params, problem_str,
                          analyzers, Ns, dts):
    """Analyze solutions for a range of discretizations of given problem with given scheme.

    Executes given problem with given scheme repeatedly in
    a loop over discretization parameters N and dt.
    After each run data is collected from all provided analyzers.

    Collected data is returned in a dict {field_name: {(N,dt): field_data}}.
    """
    # ... Solve problems for each discretization
    data = {}
    for problem_params.N in Ns:
        for problem_params.dt in dts:
            s = Scheme(scheme_params)
            p = Problem(problem_params)

            casedir = "results/%s/%s/N=%d/dt=%e" % (scheme_str, problem_str, N, dt)
            pp = NSPostProcessor({"casedir":casedir})

            fields = [Analyzer(params=analyzer_params) for (Analyzer,analyzer_params) in analyzers]
            pp.add_fields(fields)

            nssolver = NSSolver(p, s, pp)

            try:
                # Solve the problem! This will take some time...
                nssolver.solve()

                # Make dicts first time
                for field in fields:
                    if not field.name in data:
                        data[field.name] = {}

                # Collect data from analyzers
                for field in fields:
                    # TODO: What kind of interface is this get_data()? Is it ppfields or specific for testing?
                    data[field.name][(problem_params.N, problem_params.dt)] = field.get_data()

            except Exception as e:
                print "The scheme did not work (%s, %s)" % (scheme_str, problem_str)
                print e
    return data

class TestDiscretizationVariations(unittest.TestCase):

    def _analyze_spatial_convergence(self, data, Ns, dts):
        velocity = data["Velocity"]
        print velocity.keys()

        # With time discretization fixed, loop over spatial discretization and compute norms
        norms_h = {}
        for dt in dts:
            norms_h[dt] = [None]*(len(Ns)-1)

            for i in range(len(Ns)-1):
                N0 = Ns[i]
                N1 = Ns[i+1]

                try:
                    u_fine = velocity[(N1,dt)]["data"]
                    u_coarse = velocity[(N0, dt)]["data"]

                    uc = dolfin.interpolate(u_coarse, u_fine.function_space())
                    u_diff_norm = dolfin.assemble((u_fine - uc)**2*dolfin.dx())
                    u_norm = dolfin.assemble(u_fine**2*dolfin.dx())

                        norms_h[dt][i] = { "u_diff_norm": u_diff_norm, "u_norm": u_norm }

                        print "difference between level ", N1, " and ", N0, " with dt ", dt, " is ", u_diff_norm, " versus u_norm ", u_norm
                    except Exception as e:
                        print "Not able to compare", (N0, N1), dt
                        print e

        # FIXME: Add assertions to validate convergence rate from norms_h[dt][:]
        for dt in dts:
            norms = norms_h[dt]
            rates = []
            for i in range(len(Ns)-1):
                try:
                    rate = norms[i+1]["u_diff_norm"] / norms[i]["u_diff_norm"]
                except:
                    rate = None
                rates.append(rate)
            print "Convergence rates w.r.t N at dt=%g: %s" % (dt, rates)

    def _analyze_temporal_convergence(self, data, Ns, dts):
        velocity = data["Velocity"]
        print velocity.keys()

        # With spatial discretization fixed, loop over time discretization and compute norms
        norms_dt = {}
        for N in Ns:
            norms_dt[N] = [None]*(len(dts)-1)

            for i in range(len(dts)-1):
                dt0 = dts[i]
                dt1 = dts[i+1]

                try:
                    u_fine = velocity[(N,dt1)]["data"]
                    u_coarse = velocity[(N, dt0)]["data"]

                    uc = dolfin.interpolate(u_coarse, u_fine.function_space())
                    u_diff_norm = dolfin.assemble((u_fine - uc)**2*dolfin.dx())
                    u_norm = dolfin.assemble(u_fine**2*dolfin.dx())

                    norms_dt[N][i] = { "u_diff_norm": u_diff_norm, "u_norm": u_norm }

                    print "difference between dt ", dt1 , " and ", dt0, " at level ", N, " is ", u_diff_norm, " versus u_norm ", u_norm

                except Exception as e:
                    print "Not able to compare ", N, (dt0, dt1)
                    print e

        # FIXME: Add assertions to validate convergence rate from norms_dt[N][:]
        for N in Ns
            norms = norms_dt[N]
            rates = []
            for i in range(len(dts)-1):
                try:
                    rate = norms[i+1]["u_diff_norm"] / norms[i]["u_diff_norm"]
                except:
                    rate = None
                rates.append(rate)
            print "Convergence rates w.r.t dt  at N=%g: %s" % (N, rates)


    #@unittest.skip("Enable validation tests when set up properly!")
    def test_Q(self):
        # FIXME: Implement test_Q here, using the QDeltaLambda2 analyzer, or add that to test_grid_convergence.
        #        (It's not clear to me what this is, so need to check with Kent.)
        pass

    #@unittest.skip("Enable validation tests when set up properly!")
    def test_grid_convergence(self):

        # ... Choice of schemes to test
        def schemes():
            # Selecting all schemes with default parameters:
            #sch = [(Scheme, None) for Scheme in all_schemes]

            # Selecting particular schemes with default parameters:
            #sch = [(Scheme, None) for Scheme in (IPCS, IPCS_Stable, IPCS_Stabilized, SegregatedIPCS)]

            # Selecting particular schemes with custom parameters:
            sch = [
                (IPCS_Stabilized, {"theta":1.0}),
                (IPCS_Stabilized, {"theta":0.5}),
                (IPCS_Stable, {})
                ]

            for Scheme, scheme_params in sch:
                scheme_str = Scheme.shortname()
                for k in sorted(scheme_params.keys()):
                    scheme_str += "_%s=%s" % (k, str(scheme_params[k]))

                yield (Scheme, scheme_params, scheme_str)

        # ... Problem variations
        def problems():
            from flow_around_cylinder import FlowAroundACylinder
            problems = [FlowAroundACylinder]

            mus = [0.5e-2]

            for Problem in problems:
                params = Problem.default_params()
                for params.mu in mus:
                    problem_str = os.path.join(Problem.shortname(), "mu=%s" % str(params.mu))

                    yield (Problem, params, problem_str)

        # ... Configuration of postprocessing
        ppfield_pd = ParamDict(
            saveparams=ParamDict(
                save=True, # FIXME: Make storing configurable, better in an automated test to have automatic in-memory analysis
                ),
            timeparams=ParamDict(
                step_frequency=10,
                )
            )
        analyzers = [(Velocity, ppfield_pd)]
        # TODO: Generalize convergence test on analyzer classes:
        #analyzers = [(Velocity, ppfield_pd), (QDeltaLambda2, ppfield_pd)]

        # ... Discretization refinements
        #Ns = [2, 4]
        Ns = [64, 128, 256]
        #dts = [0.1, 0.05]
        dts = [0.025, 0.0125]

        # ... Analyze convergence for each scheme/problem combination
        for Scheme, scheme_params, scheme_str in schemes():
            for Problem, problem_params, problem_str in problems():
               data = run_convergence_sweep(Scheme, scheme_params, scheme_str,
                                            Problem, problem_params, problem_str,
                                            analyzers, Ns, dts)
               print ""
               print "scheme ", scheme_str
               print "problem ", problem_str
               self._analyze_spatial_convergence(data, Ns, dts)
               self._analyze_temporal_convergence(data, Ns, dts)


    def _print_table(self, data):
        # TODO: Split into computing and presenting table
        fieldname = "AnalyticalSolutionAnalyzer"
        subfields = ["u0", "u1", "u2", "p"]
        for x in subfields:
            print ""
            print ""
            print x
            print "dt ",
            for dt in dts:
                print dt,
            for N in Ns:
                print "\nN ", N,
                for dt in dts:
                    print " %2.2e " % data[fieldname][(N, dt)]["data"][x],


    #@unittest.skip("Enable validation tests when set up properly!")
    def test_analytical_solution(self):

        # ... Choice of schemes to test
        def schemes():
            # Selecting all schemes with default parameters:
            #sch = [(Scheme, None) for Scheme in all_schemes]

            # Selecting particular schemes with default parameters:
            #sch = [(Scheme, None) for Scheme in (IPCS, IPCS_Stable, IPCS_Stabilized, SegregatedIPCS)]

            # Selecting particular schemes with custom parameters:
            sch = [
                (IPCS_Stable, {})
                ]

            for Scheme, scheme_params in sch:
                scheme_str = Scheme.shortname()
                for k in sorted(scheme_params.keys()):
                    scheme_str += "_%s=%s" % (k, str(scheme_params[k]))

                yield (Scheme, scheme_params, scheme_str)

        # ... Problem variations
        def problems():
            from beltrami import Beltrami
            problems = [Beltrami]

            for Problem in problems:
                params = Problem.default_params()
                problem_str = Problem.shortname()

                yield (Problem, params, problem_str)

        # ... Configuration of postprocessing
        # FIXME: Make storing configurable, better in an automated test to have automatic in-memory analysis
        ppfield_pd1 = ParamDict(
            saveparams=ParamDict(
                save=True,
                ),
            )
        ppfield_pd3 = ParamDict(
            saveparams=ParamDict(
                save=True,
                ),
            timeparams=ParamDict(
                step_frequency=10,
                )
            )
        analyzers = [
            (AnalyticalSolutionAnalyzer, ppfield_pd1),
            (EnergyAnalyzer, ppfield_pd1),
            (WSS, ppfield_pd3),
            ]

        # ... Discretization refinements
        Ns = [2, 4, 8, 16]
        dts = [0.1, 0.05, 0.025, 0.0125]

        # ... Analyze convergence for each scheme/problem combination
        for Scheme, scheme_params, scheme_str in schemes():
            for Problem, problem_params, problem_str in problems():
               data = run_convergence_sweep(Scheme, scheme_params, scheme_str,
                                            Problem, problem_params, problem_str,
                                            analyzers, Ns, dts)
               print ""
               print "scheme ", scheme_str
               print "problem ", problem_str
               self._print_table(data)


    #@unittest.skip("Enable validation tests when set up properly!")
    def test_analytical_solution_stabilized(self):

        # ... Choice of schemes to test
        def schemes():
            # Selecting all schemes with default parameters:
            #sch = [(Scheme, None) for Scheme in all_schemes]

            # Selecting particular schemes with default parameters:
            #sch = [(Scheme, None) for Scheme in (IPCS, IPCS_Stable, IPCS_Stabilized, SegregatedIPCS)]

            # Selecting particular schemes with custom parameters:
            sch = [
                (IPCS_Stabilized, {"theta":1.0}),
                (IPCS_Stabilized, {"theta":0.5}),
                (IPCS_Stable, {})
                ]

            for Scheme, scheme_params in sch:
                scheme_str = Scheme.shortname()
                for k in sorted(scheme_params.keys()):
                    scheme_str += "_%s=%s" % (k, str(scheme_params[k]))

                yield (Scheme, scheme_params, scheme_str)

        # ... Problem variations
        def problems():
            from beltrami import Beltrami
            problems = [Beltrami]

            mus = [1.0, 1.0e-3]
            mus = [1.0e-3]

            for Problem in problems:
                params = Problem.default_params()
                for params.mu in mus:
                    problem_str = os.path.join(Problem.shortname(), "mu=%s" % str(params.mu))

                    yield (Problem, params, problem_str)

        # ... Configuration of postprocessing
        # FIXME: Make storing configurable, better in an automated test to have automatic in-memory analysis
        ppfield_pd1 = ParamDict(
            saveparams=ParamDict(
                save=True,
                ),
            )
        ppfield_pd3 = ParamDict(
            saveparams=ParamDict(
                save=True,
                ),
            timeparams=ParamDict(
                step_frequency=10,
                )
            )
        analyzers = [
            (AnalyticalSolutionAnalyzer, ppfield_pd1),
            (EnergyAnalyzer, ppfield_pd1),
            (WSS, ppfield_pd3),
            ]

        # ... Discretization refinements
        #Ns = [2, 4]
        Ns = [2, 4, 8]
        #dts = [0.1, 0.05]
        #dts = [0.5]
        dts = [0.1, 0.05]

        # ... Analyze convergence for each scheme/problem combination
        for Scheme, scheme_params, scheme_str in schemes():
            for Problem, problem_params, problem_str in problems():
               data = run_convergence_sweep(Scheme, scheme_params, scheme_str,
                                            Problem, problem_params, problem_str,
                                            analyzers, Ns, dts)
               print ""
               print "scheme ", scheme_str
               print "problem ", problem_str
               self._print_table(data)
