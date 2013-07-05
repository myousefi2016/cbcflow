import sys
sys.path.insert(0, "../../site-packages/")


from headflow import *
from headflow.dol import *

from headflow.core.spaces import NSSpacePool

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
            #exit()



mesh = Mesh("dog_mesh_37k.xml.gz")

problem = NSProblem()
problem.mesh = mesh
problem.params.T = 10

V = FunctionSpace(mesh, "CG", 1)

slicemesh = create_slicemesh([55, 44, 29], [0,-1,0], [1,0,0], 10, 10)
boxmesh = BoxMesh(50, 39, 24, 60, 49, 34, 20,20,20)
#submesh = BoxMesh(50, 39, 24, 60, 49, 34, 10,10,10)
#submesh = BoxMesh(35, 24, 9, 75, 64, 49, 20,20,20)
subfunc1 = SubFunction("Velocity", slicemesh, {"save": True})
subfunc2 = SubFunction("Pressure", boxmesh, {"save": True})


pp = NSPostProcessor()
#pp.add_fields([subfunc1, subfunc2])
pp.add_field(subfunc1)

scheme = Scheme()

solver = NSSolver(problem, scheme, pp)
solver.solve()


#subfunction = SubFunction()