#!/usr/bin/env python


from cbcflow import *
from cbcflow.dol import *

from math import pi, e

class Beltrami(NSProblem):
    "3D test problem with known analytical solution."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Create mesh of box (-1, 1) x (-1, 1) x (-1, 1)
        refinement_levels = [4,8,16,32,64]
        N = refinement_levels[self.params.refinement_level]
        mesh = UnitCubeMesh(N, N, N)
        scaled = 2*(mesh.coordinates() - 0.5)
        mesh.coordinates()[:, :] = scaled
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(1)
        DomainBoundary().mark(facet_domains, 0)

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=0.5,
            dt=0.05,
            # Physical parameters
            rho=1.0,
            mu=1.0,
            )
        params.update(
            # Spatial parameters
            refinement_level=0,
            )
        return params

    def analytical_solution(self, spaces, t):
        # The analytical solution
        # Velocity
        analytical_u = \
            ('-((a*(pow(E,a*x[2])*cos(a*x[0] + d*x[1]) + pow(E,a*x[0])*sin(a*x[1] + d*x[2])))/pow(E,pow(d,2)*t*nu))',
             '-((a*(pow(E,a*x[0])*cos(a*x[1] + d*x[2]) + pow(E,a*x[1])*sin(d*x[0] + a*x[2])))/pow(E,pow(d,2)*t*nu))',
             '-((a*(pow(E,a*x[1])*cos(d*x[0] + a*x[2]) + pow(E,a*x[2])*sin(a*x[0] + d*x[1])))/pow(E,pow(d,2)*t*nu))')

        # Pressure
        analytical_p = \
            ('-(1.0/2.0)*(pow(a,2)*(pow(E,2*a*x[0]) + pow(E,2*a*x[1]) + pow(E,2*a*x[2]) + 2*pow(E,a*(x[1] + x[2]))*cos(d*x[0] + a*x[2])*sin(a*x[0] + d*x[1]) + 2*pow(E,a*(x[0] + x[1]))*cos(a*x[1] + d*x[2])*sin(d*x[0] + a*x[2]) + 2*pow(E,a*(x[0] + x[2]))*cos(a*x[0] + d*x[1])*sin(a*x[1] + d*x[2])))/(pow(E,pow(d,2)*t*nu))')

        # Common parameters pertinent to the functional forms above
        u_params = {'a': pi/4.0, 'd': pi/2.0, 'E': e,                         'nu': self.params.mu/self.params.rho, 't': t}
        p_params = {'a': pi/4.0, 'd': pi/2.0, 'E': e, 'rho': self.params.rho, 'nu': self.params.mu/self.params.rho, 't': t}

        # Compile expressions
        exact_u = [Expression(analytical_u[d], **u_params) for d in xrange(3)]
        exact_p = Expression(analytical_p, **p_params)
        
        return [as_vector(exact_u), exact_p]
    
    def test_references(self, spaces, t):
        return self.analytical_solution(spaces, t)
    
    def test_fields(self):
        return [SolutionField("Velocity"), SolutionField("Pressure")]
    

    def initial_conditions(self, spaces, controls):
        #exact_u, exact_p = self.analytical_solution(spaces, t=0.0)
        exact_u, exact_p = self.test_references(spaces, 0.0)
        self.p = Function(spaces.Q)
        return (exact_u, exact_p)

    def boundary_conditions(self, spaces, u, p, t, controls):
        exact_u, exact_p = self.test_references(spaces, float(t))
        bcu = [(exact_u, 0)]
        bcp = []
        return bcu, bcp

    def update(self, spaces, u, p, t, timestep, bcs, observations, controls):
        bcu, bcp = bcs
        uve = bcu[0][0]
        for ue in uve: ue.t = float(t)
        

    '''
    # FIXME: Change this to use the new test_functionals, test_references interface:
    def functional(self, t, u, p):
        if t < self.T:
            return 0.0
        else:
            exact_u, exact_p = self.analytical_solution(spaces, t=t)

            error = 0
            for exact_u, calc_u in zip(exact_u, u):
                error += sqr(errornorm(exact_u, calc_u) / norm(exact_u, mesh=self.mesh))

            return sqrt(error / len(u))

    def reference(self, t):
        return 0.0
    '''

def main():
    set_log_level(100)
    problem = Beltrami()
    scheme = IPCS_Stable()

    casedir = "results_demo_%s_%s" % (problem.shortname(), scheme.shortname())
    plot_and_save = dict(plot=True, save=True)
    fields = [
        SolutionField("Pressure", plot_and_save),
        SolutionField("Velocity", plot_and_save),
        ]
    postproc = PostProcessor({"casedir": casedir})
    postproc.add_fields(fields)

    solver = NSSolver(problem, scheme, postproc)
    solver.solve()

if __name__ == "__main__":
    main()
