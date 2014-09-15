.. _Womersley3D:

Womersley flow in 3D
===================================

In this demo it is demonstrated how to handle problems with time-dependent boundary
conditions and known analytical solution/reference solution. The problem is transient
Womersley flow in a cylindrical pipe.

The source code can be found in :download:`Womersley3D.py`.

We start by importing cbcflow and dolfin: ::


    from cbcflow import *
    from dolfin import *
    

Specifying the domain
____________________________________

Our domain is a cylinder of length 10.0 and radius 0.5: ::

    LENGTH = 10.0
    RADIUS = 0.5
    
The meshes for this has been pregenerated and is available in the demo data, see :ref:`Demos`. ::
    
    files = [path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/pipe_1k.xml.gz"),
             path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/pipe_3k.xml.gz"),
             path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/pipe_24k.xml.gz"),
             path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/pipe_203k.xml.gz"),
             path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/pipe_1611k.xml.gz"),
            ]
    
We define *SubDomain* classes for inflow and outflow: ::
    
    class Inflow(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] < 1e-6 and on_boundary
    
    class Outflow(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] > LENGTH-1e-6 and on_boundary

We could also add a *SubDomain*-class for the remaining wall, but this will be handled later.

Defining the NSProblem
_____________________________________

We first define a problem class inheriting from :class:`.NSProblem`: ::
    
    class Womersley3D(NSProblem):
    
The parameters of the problem are defined to give a Reynolds number of about 30 and
Womersley number of about 60. ::

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=None,
            dt=1e-3,
            period=0.8,
            num_periods=1.0,
            # Physical parameters
            rho=1.0,
            mu=1.0/30.0,
            )
        params.update(
            # Spatial parameters
            refinement_level=0,
            # Analytical solution parameters
            Q=1.0,
            )
        return params

In the constructor, we load the mesh from file and mark the boundary domains relating
to inflow, outflow and wall in a *FacetFunction*: ::

    def __init__(self, params=None):
        NSProblem.__init__(self, params)
        
        # Load mesh
        mesh = Mesh(files[self.params.refinement_level])

        # We know that the mesh contains markers with these id values
        self.wall_boundary_id = 0
        self.left_boundary_id = 1
        self.right_boundary_id = 2
        
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(3)
        DomainBoundary().mark(facet_domains, self.wall_boundary_id)
        Inflow().mark(facet_domains, self.left_boundary_id)
        Outflow().mark(facet_domains, self.right_boundary_id)
        
We then define a transient profile for the flow rate, for use later: ::

         # Setup analytical solution constants
        Q = self.params.Q
        self.nu = self.params.mu / self.params.rho

        # Beta is the Poiseuille pressure drop if the flow rate is stationary Q
        self.beta = 4.0 * self.nu * Q / (pi * RADIUS**4)

        # Setup transient flow rate coefficients
        print "Using transient bcs."
        P = self.params.period
        tvalues = np.linspace(0.0, P)
        Qfloor, Qpeak = -0.2, 1.0
        Qvalues = Q * (Qfloor + (Qpeak-Qfloor)*np.sin(pi*((P-tvalues)/P)**2)**2)
        self.Q_coeffs = zip(tvalues, Qvalues)
        
Finally, we store the mesh and facet domains to *self*: ::

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

The analytical solution
___________________________________

The Womersley profile can be obtained by using the helper function :func:`.make_womersley_bcs`.
This function returns a list of scalar *Expression* instances defining the Womersley profile: ::

    def analytical_solution(self, spaces, t):
        # Create womersley objects
        ua = make_womersley_bcs(self.Q_coeffs, self.mesh, self.left_boundary_id, self.nu, None, self.facet_domains)
        for uc in ua:
            uc.set_t(t)
        pa = Expression("-beta * x[0]", beta=1.0)
        pa.beta = self.beta # TODO: This is not correct unless stationary...
        return (ua, pa)
        
Note that the pressure solution defined here is not correct in the transient case.

Using an analytical/reference solution
_________________________________________

If one for example wants to validate a scheme, it is required to define the following functions: ::
    
    def test_fields(self):
        return [Velocity(), Pressure()]
    
    def test_references(self, spaces, t):
        return self.analytical_solution(spaces, t)
        
The :func:`test_fields` function tells that the fields :class:`.Velocity` and
:class:`.Pressure` should be compared to the results from :func:`test_references`, namely
the analytical solution.

These functions are used in the regression/validation test suite to check and record errors.

Initial conditions
____________________________________

As initial conditions we simply use the analytical solution at t=0.0: ::

    def initial_conditions(self, spaces, controls):
        return self.analytical_solution(spaces, 0.0)

Boundary conditions
____________________________________

At the boundaries, we also take advantage of the analytical solution, and we set no-slip
conditions at the cylinder walls:

    def boundary_conditions(self, spaces, u, p, t, controls):
        # Create no-slip bcs
        d = len(u)
        u0 = [Constant(0.0)] * d
        noslip = (u0, self.wall_boundary_id)

        # Get other bcs from analytical solution functions
        ua, pa = self.analytical_solution(spaces, t)

        # Create inflow boundary conditions for velocity
        inflow = (ua, self.left_boundary_id)

        # Create outflow boundary conditions for pressure
        p_outflow = (pa, self.right_boundary_id)

        # Return bcs in two lists
        bcu = [noslip, inflow]
        bcp = [p_outflow]
        
        return (bcu, bcp)
        
Now, since these boundary conditions are transient, we need to use the :func:`update` function.
The :func:`boundary_conditions` function is called at the start of solve step, and a
call-back is done to the :func:`update` function to do any updates to for example the
boundary conditions. In here, we update the time in the inlet boundary condition: ::

    def update(self, spaces, u, p, t, timestep, bcs, observations, controls):
        bcu, bcp = bcs
        uin = bcu[1][0]
        for ucomp in uin:
            ucomp.set_t(t)

Solving the problem
________________________________

Finally, we initate the problem, a scheme and postprocessor ::
    
    def main():
        problem = Womersley3D({"refinement_level": 2})
        scheme = IPCS_Stable()
    
        casedir = "results_demo_%s_%s" % (problem.shortname(), scheme.shortname())
        plot_and_save = dict(plot=True, save=True)
        fields = [
            Pressure(plot_and_save),
            Velocity(plot_and_save),
            ]
        postproc = PostProcessor({"casedir": casedir})
        postproc.add_fields(fields)

and solves the problem ::
    
        solver = NSSolver(problem, scheme, postproc)
        solver.solve()
    