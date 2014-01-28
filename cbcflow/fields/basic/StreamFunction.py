from ..bases.PPField import PPField
from dolfin import TrialFunction, TestFunction, dot, grad, DirichletBC, DomainBoundary, dx, Constant, assemble, KrylovSolver, Vector, Function

class StreamFunction(PPField):
    def before_first_compute(self, pp, spaces, problem):
        u = pp.get("Velocity")
        assert len(u) == 2, "Can only compute stream function for 2D problems"
        V = spaces.DU0
        psi = TrialFunction(V)
        self.q = TestFunction(V)
        a = dot(grad(psi), grad(self.q))*dx()
        
        
        self.bc = DirichletBC(V, Constant(0), DomainBoundary())
        A = assemble(a)
        self.L = Vector()
        self.bc.apply(A)
        self.solver = KrylovSolver(A, "cg")
        self.psi = Function(V)
        

    def compute(self, pp, spaces, problem):
        u = pp.get("Velocity")
        assemble(dot(u[1].dx(0)-u[0].dx(1), self.q)*dx(), tensor=self.L)
        self.bc.apply(self.L)
        
        self.solver.solve(self.psi.vector(), self.L)
        
        return self.psi
        
        
        
        
        
        
     
