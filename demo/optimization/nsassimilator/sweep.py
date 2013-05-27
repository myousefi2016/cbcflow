
from run import default_params, run

# ====== Iterator over all parameter combinations

def sweep():
    verbose = True
    k = 0
    casedir = "testcase%d"
    #for dt in (0.1, 0.05, 0.025, 0.0125):
    for dt in (0.05,):
        #for mu in (1.0, 0.035):
        for mu in (1.0,):
            #for alpha in (1e-3, 0.0):
            for alpha in (0.0,):
                #for scale in ("auto", 1e4):
                for scale in ("auto",):
                    #for pdim in (1, 5, 9, 17, 33):
                    for pdim in (9,):
                        #for minimize_tolerance in (1e-4, 1e-6):
                        for minimize_tolerance in (1e-6,):

                            p = default_params()
                            pp = p.problem
                            sp = p.scheme
                            ap = p.assimilator

                            #sp.form_compiler_parameters={'quadrature_degree': 3}
                            sp.verbose = verbose

                            pp.alpha = alpha
                            pp.dt = dt
                            pp.pdim = pdim
                            pp.scale = scale
                            pp.mu = mu
                            pp.num_periods = 1

                            pp.num_timesteps = 1 # TEMPORARY LIMIT FOR TESTING

                            ap.minimize_tolerance = minimize_tolerance
                            ap.verbose = verbose
                            ap.casedir = casedir % k

                            k += 1
                            yield p


# ====== Problem independent script stuff

def pick_params(i):
    p = None
    for j, q in enumerate(sweep()):
        if i == j:
            p = q
            break
    return p

def show_parameter_sets():
    k = 0
    for p in sweep():
        k += 1
        print "------ Parameter set %d:" % k
        print p
    print "------ Total %d parameter sets" % k

def main(args):
    a = args[0]

    if len(args) == 1:
        p = None
    else:
        i = int(args[1])
        p = pick_params(i)
        if p is None:
            print "Invalid parameter set number %d" % i
            return 0

    if a == "show":
        if p is None:
            show_parameter_sets()
        else:
            print "Showing parameter set %d:" % i
            print p
            return 0
        return 0
    elif a == "run":
        if p is None:
            for i, p in enumerate(sweep()):
                print "Running with parameter set %d:" % i
                print p
                run(p)
        else:
            print "Running with parameter set %d:" % i
            print p
            run(p)
        return 0
    else:
        print "Unknown action %s, try 'show' or 'run'." % a
        return 1

if __name__ == "__main__":
    import sys
    args = sys.argv[1:]
    sys.exit(main(args))
