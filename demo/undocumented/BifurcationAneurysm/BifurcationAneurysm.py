#!/usr/bin/env python

from cbcflow import *
from cbcflow.dol import *
from os import path

files = [path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/dog_mesh_37k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"../../../cbcflow-data/dog_mesh_97k.xml.gz"),
        ]

c0 = Constant(0)

class BifurcationAneurysm(NSProblem):
    "Template of a typical implementation of a NSProblem subclass"

    def __init__(self, params=None):
        """Initialize problem. The following are required:
        - self.mesh
        etc.
        """
        NSProblem.__init__(self, params)

        mesh = Mesh(files[self.params.refinement_level])
        
        cell_domains = MeshFunction("size_t", mesh, 3)
        cell_domains.set_all(0)
        subdomain = CompiledSubDomain("x[0] > 49.0 && x[0] < 60.0 \
                                      && x[1] > 37.4 && x[1] < 48.9 \
                                      && x[2] > 25.2 && x[2] < 40.4")

        subdomain.mark(cell_domains, 1)
        
        self.initialize_geometry(mesh, cell_domains=cell_domains)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            dt=1e-2,
            #period=1.0,
            #num_periods=0.3,
            T=2,
            # Physical parameters
            mu=0.00345,
            rho=0.00106,
            )
        params.update(
            # Spatial discretization parameters
            refinement_level=0,
            )
        return params

    def initial_conditions(self, spaces, controls):
        "Return initial conditions as list of scalars (velocity) and scalar (pressure)."
        icu = [c0, c0, c0]
        icp = c0
        return (icu, icp)

    def boundary_conditions(self, spaces, u, p, t, controls):
        "Return boundary conditions as lists."

        factor = 1000
        profile = [0.4, 1.6, 1.4, 1.0, 0.8, 0.6, 0.55, 0.5, 0.5, 0.45, 0.4]
        profile = [p*factor for p in profile]
        time = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

        #inflow = Womersley(zip(time, profile), self.mesh, 1, self.mu/self.rho, scale_to=1000)
        #inflow = Poiseuille(zip(time, profile), self.mesh, 1, scale_to=1000)

        inflow = Poiseuille(zip(time, profile), self.mesh, 1)
        for e in inflow: e.set_t(float(t))

        bcu = [(inflow, 1),
               ([c0, c0, c0], 0)]
        bcp = [(c0, 2),
               (c0, 3)]
        return (bcu, bcp)

    def update(self, spaces, u, p, t, timestep, bcs, observations, controls):
        bcu, bcp = bcs
        inflow = bcu[0][0]
        for e in inflow: e.set_t(float(t))
        
        


def main():
    set_log_level(60)
    dt = 1e-2
    problem = BifurcationAneurysm(dict(refinement_level=0, dt=dt, T=2.0))
    print problem.mesh
    scheme = IPCS_Stable(dict(
        rebuild_prec_frequency = 1,
        u_tent_prec_structure = "same_nonzero_pattern",
        #p_corr_solver_parameters = dict(relative_tolerance=1e-6, absolute_tolerance=1e-6, monitor_convergence=False),
        u_degree=1,
        p_degree=1,
        solver_u_tent=("gmres", "additive_schwarz"),
        #solver_u_corr=("cg", "additive_schwarz"),
        solver_u_corr = "WeightedGradient",
        theta=1.0,
        ))
    
    parameters["krylov_solver"]["relative_tolerance"] = 1e-15
    parameters["krylov_solver"]["absolute_tolerance"] = 1e-15

    casedir = "results_demo_%s_%s" % (problem.shortname(), scheme.shortname())
    plot_and_save = dict(plot=False, save=True, stride_timestep=100)
    
    postproc = NSPostProcessor({"casedir": "Results/"})
    
    
    def set_up_fields(problem):
        T0 = problem.params.T/2.0
        T1 = problem.params.T
        fields = []
        
        plot = True
        
        # Basic fields
        fields.append(Pressure(dict(plot=plot, save=True, stride_timestep=10)))
        fields.append(Velocity(dict(plot=plot, save=True, stride_timestep=10)))
    
        # On boundary
        fields.append(WSS(dict(plot=plot, save=True, start_time=0)))
        fields.append(Boundary("Pressure", dict(plot=plot, save=True)))
        
        # Time-integrated fields
        fields.append(TimeIntegral("WSS", dict(save=True, start_time=T0, end_time=T1)))
        fields.append(OSI(dict(save=True, start_time=T0, end_time=T1)))
        
        # SubFunctions and Restrictions
        from cbcflow.utils.common.submesh import create_submesh
        from cbcflow.utils.fields.Slice import Slice
        
        submesh = create_submesh(problem.mesh, problem.cell_domains, 1)
        fields.append(Restrict("Velocity", submesh, dict(save=True, plot=plot)))

        # Need mpi4py for SubFunction        
        try:
            import mpi4py
        except:
            mpi4py = None
        
        if mpi4py != None:          
            slicemesh = Slice(problem.mesh, (54.8, 44.0, 33.1), (-0.23, -0.10, 0.97))
            fields.append(SubFunction("Velocity", slicemesh, dict(save=True, plot=plot, plot_args=dict(mode='color'))))
            fields.append(SubFunction("Pressure", slicemesh, dict(save=True, plot=plot, plot_args=dict(mode='color'))))
            
        # Point evaluation
        fields.append(PointEval("Velocity", ((54.8, 44.0, 33.1), (53.6, 43.6, 37.9)), dict(save=True)))
        
        # Derivatives
        fields.append(TimeDerivative("Pressure", dict(save=True, plot=plot)))
        
        return fields

    postproc.add_fields(set_up_fields(problem))
    
    solver = NSSolver(problem, scheme, postproc, dict(timer_frequency=10))
    solver.solve()
    

if __name__ == "__main__":
    main()
