from .MetaPPField import MetaPPField
from dolfin import assemble, dx, Function, Constant

class DomainAvg(MetaPPField):
    def compute(self, pp, spaces, problem):
        u = pp.get(self.valuename)

        # Find mesh/domain
        if isinstance(u, Function):
            mesh = u.function_space().mesh()
        else:
            mesh = problem.mesh
        
        # Calculate volume
        if not hasattr(self, "volume"):
            self.volume = assemble(Constant(1)*dx(), mesh=mesh)
        
        if u.rank() == 0:
            value = assemble(u*dx(), mesh=mesh)/self.volume
        elif u.rank() == 1:
            value = [assemble(u[i]*dx(), mesh=mesh)/self.volume for i in xrange(u.value_size())]
        elif u.rank() == 2:
            value = []
            for i in xrange(u.shape()[0]):
                for j in xrange(u.shape()[1]):
                    value.append(assemble(u[i,j]*dx(), mesh=mesh)/self.volume)

        return value
