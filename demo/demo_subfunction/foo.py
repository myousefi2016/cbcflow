import sys
sys.path.insert(0, "../../site-packages/")


from headflow import *
from headflow.dol import *

from headflow.core.spaces import NSSpacePool
set_log_level(100)
parameters["allow_extrapolation"] = True
class Scheme(NSScheme):
    def solve(self, problem, update):
        mesh = problem.mesh
        self.params.u_degree = 1
        self.params.p_degree = 1
        spaces = NSSpacePool(mesh, self.params.u_degree, self.params.p_degree)
        Q = spaces.Q
        V = spaces.V
        for i in range(1,5):
            p = Function(Q, "pressure%d.xml.gz" %i)
            u = Function(V, "velocity%d.xml.gz" %i)
            update(u, p, i*0.1, i, spaces)
            
        return {"spaces": spaces}



mesh = Mesh("dog_mesh_37k.xml.gz")

problem = NSProblem({"T": 10.0})
problem.mesh = mesh
problem.params.mu = 1.0
    
V = FunctionSpace(mesh, "CG", 1)

# Create slice (basemesh, origin, normal)
slicemesh1 = Slice(problem.mesh, [55, 44, 29], [0,0,1])
slicemesh2 = Slice(problem.mesh, [50, 44, 20], [0,0,1])
slicemesh3 = Slice(problem.mesh, [60, 44, 25], [0,0,1])
boxmesh = BoxMesh(50, 39, 24, 60, 49, 34, 20,20,20)

u = Velocity({"save": True})
p = Pressure({"save": True})

subfunc1 = SubFunction("Velocity", slicemesh1, {"save": True}, label=1)
subfunc2 = SubFunction("Velocity", slicemesh2, {"save": True}, label=2)
subfunc3 = SubFunction("Pressure", slicemesh3, {"save": True}, label=3)
subfunc4 = SubFunction("Pressure", boxmesh, {"save" : True})

davg1 = DomainAvg(subfunc1, {"save": True})
davg2 = DomainAvg(subfunc2, {"save": True})
davg3 = DomainAvg(subfunc3, {"save": True})
davg4 = DomainAvg(subfunc4, {"save": True})


pp = NSPostProcessor()
pp.add_fields([subfunc1, subfunc2, subfunc3, subfunc4, davg1, davg2, davg3, davg4])

scheme = Scheme()

solver = NSSolver(problem, scheme, pp)
solver.solve()
