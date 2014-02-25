We start by importing cbcflow and dolfin: ::

   from cbcflow import *
   from dolfin import *

Specifying the domain
____________________________________

The meshes for this problem is pregenerated, and is specified at the following locations: ::

    from os import path

    files = [path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/cylinder_0.6k.xml.gz"),
             path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/cylinder_2k.xml.gz"),
             path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/cylinder_8k.xml.gz"),
             path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/cylinder_32k.xml.gz"),
             path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/cylinder_129k.xml.gz"),
            ]

This requires that you have installed the demo data, as specified in :ref:`getstarted`.

The domain is based on a rectangle with corners in (0,0), (0,1), (10,0) and (10,1).
The cylinder is centered in (2,0.5) with radius of 0.12. The different boundaries
of the domain is specified as: ::

    class LeftBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], 0.0)
    
    class RightBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], 10.0)
    
    class Cylinder(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and (sqrt((x[0]-2.0)**2+(x[1]-0.5)**2) < 0.12+DOLFIN_EPS)
    
    class Wall(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and (near(x[1], 0.0) or near(x[1], 1.0))

Defining a NSProblem
__________________________________

To define a problem class recognized by cbcflow, the class must inherit from
:class:`.NSProblem`:

.. code-block:: python

    class FlowAroundCylinder(NSProblem):
    
Parameters
--------------------------------------
This class inherit from the :class:`.Parameterized` class,
allowing for parameters in the class interface. We supply default parameters to
the problem: ::

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=5.0,
            dt=0.1,
            # Physical parameters
            rho=1.0,
            mu=1.0/1000.0,
            )
        params.update(
            # Spatial parameters
            refinement_level=0,
            )
        return params
        
This takes the default parameters from NSProblem and replaces some parameters common
for all NSProblems. We set the end time to 5.0 with a timestep of 0.1, the density
:math:`\rho=1.0` and dynamic viscosity :math:`\mu=0.001`. In addition, we add
a new parameter, refinement_level, to determine which of the previously specified
mesh files to use.

Constructor
-----------------------------------
To initiate a FlowAroundCylinder-instance, we load the mesh and initialize the
geometry: ::

    def __init__(self, params=None):
        NSProblem.__init__(self, params)
        
        # Load mesh
        mesh = Mesh(files[self.params.refinement_level])

        # Create boundary markers
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(4)
        Wall().mark(facet_domains, 0)
        Cylinder().mark(facet_domains, 0)
        LeftBoundary().mark(facet_domains, 1)
        RightBoundary().mark(facet_domains, 2)

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)
        
The first call to NSProblem.__init__ updates the default parameters with any parameters
passed to the constructor as a dict or
:class:`.ParamDict`. This sets params as an
attribute to self. We load the mesh from a string defined in the files-list, and define
its domains. Finally, we call self.initialize_geometry to attach facet_domains to the mesh,
and the mesh to self.

Initial conditions
-----------------------------------
At the initial time, the fluid is set to rest, with a zero pressure gradient.
These initial conditions are prescribed by ::

    def initial_conditions(self, spaces, controls):
        c0 = Constant(0)
        u0 = [c0, c0]
        p0 = c0
        return (u0, p0)
        
The argument *spaces* is a :class:`.NSSpacePool`
helper object used to construct and contain the common function spaces related
to the Navier-Stokes solution. This is used to limit the memory consumption and
simplify the interface, so that you can, for example, call spaces.DV to get the
tensor valued gradient space of the velocity regardless of velocity degree.

The argument *controls* is used for adjoint problems, and can be disregarded for
simple forward problems such as this.


Boundary conditions
------------------------------------
As boundary conditions, we set no-slip conditions on the cylinder, at y=0.0 and y=1.0.
At the inlet we set a uniform velocity of (1.0,0.0), and zero-pressure boundary
condition at the outlet.

To determine domain to apply boundary condition, we utilize the definition of
*facet_domains* from the constructor. ::

    def boundary_conditions(self, spaces, u, p, t, controls):
        c0 = Constant(0)
        c1 = Constant(1)

        # Create no-slip boundary condition for velocity
        bcu0 = ([c0, c0], 0)
        bcu1 = ([c1, c0], 1)

        # Create boundary conditions for pressure
        bcp0 = (c0, 2)

        # Collect and return
        bcu = [bcu1, bcu2]
        bcp = [bcp0]
        return (bcu, bcp)

The way these boundary conditions are applied to the equations are determined by
the scheme used to solve the equation.

Setting up the solver
_____________________________________

Now that our *FlowAroundCylinder*-class is sufficiently defined, we can start
thinking about solving our equations. We start by creating an instance of 
*FlowAroundCylinder* class: ::

    problem = FlowAroundCylinder({"refinement_level": 2})

Note that we can pass a dict to the constructor to set, in this example, the
desired refinement level of our mesh.
    
Selecting a scheme
--------------------------------------
Several schemes are implemented in cbcflow, but only a couple are properly tested
and validated, and hence classified as *official*. Use ::

    show_schemes()

to list all schemes available, both official and unofficial.

In our application we select a very efficient operator-splitting scheme, :class:`.IPCS_Stable`, ::

    scheme = IPCS_Stable()

Setting up postprocessing
--------------------------------------
The postprocessing is set up to determine what we want to do with our obtained solution.
We start by creating a
:class:`.NSPostProcessor`
to handle all the logic: ::

    casedir = "results_demo_%s_%s" % (problem.shortname(), scheme.shortname())
    postprocessor = NSPostProcessor({"casedir": casedir})

The *casedir* parameter points the postprocessor to the directory where it should save
the data it is being asked to save. By default, it stores the mesh, all parameters and
a *play log* in that directory.

Then, we have to choose what we want to compute from the solution. The command ::

    show_fields()

lists all available :class:`.PPField`
to compute from the solution.

In this case, we are interested in the velocity, pressure and stream function,
and we wish to both plot and save these at every timestep: ::

    plot_and_save = dict(plot=True, save=True)
    fields = [
        Pressure(plot_and_save),
        Velocity(plot_and_save),
        StreamFunction(plot_and_save),
        ]

With no saveformat prescribed, the postprocessor will choose default saveformats based
on the type of data. You can use ::

    print PPField.default_parameters()

to see common parameters of these fields.

Finally, we need to add these fields to the postprocessor: ::
   
    postprocessor.add_fields(fields)
    
    
Solving the problem
----------------------------------------
We now have instances of the classes
:class:`.NSProblem`,
:class:`.NSScheme`,
and :class:`.NSPostProcessor`.

These can be combined in a general class to handle the logic between the classes,
namely a :class:`.NSSolver` instance: ::

    solver = NSSolver(problem, scheme, postprocessor)

This class has functionality to pass the solution from scheme on to the postprocessor,
report progress to screen and so on. To solve the problem, simply execute ::

    solver.solve()

    

