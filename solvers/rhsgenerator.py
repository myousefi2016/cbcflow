from dolfin import *
from ufl.form import Form

class RhsGenerator(object):
    """The instructions to create b."""
    def __init__(self, space):
        self.space = space
        self.matvecs = []
        self.forms = []
        self.vecs = []

    def __iadd__(self, ins):
        if isinstance(ins, tuple):
            A, x = ins
            if isinstance(A, Form):
                A = assemble(A)
            self.matvecs.append((A, self._as_vector_or_timedep(x), 1))
        elif isinstance(ins, GenericVector):
            self.vecs.append(ins)
        elif isinstance(ins, Form):
            self.forms.append(ins)
        else:
            raise RuntimeError, "Unknown RHS generator "+str(type(ins))
        return self

    def __isub__(self, ins):
        if isinstance(ins, tuple):
            A, x = ins
            if isinstance(A, GenericMatrix):
                self.matvecs.append((A, self._as_vector_or_timedep(x), -1))
                return self
        raise RuntimeError, "Try '+=' instead"

    def _as_vector_or_timedep(self, x):
        if isinstance(x, (GenericVector, Expression, Function)):
            return x
        return assemble(inner(x, TestFunction(self.space)) * dx)

    def _as_vector(self, x):
        if isinstance(x, GenericVector):
            return x
        if isinstance(x, Function):
            return x.vector()
        return assemble(inner(x, TestFunction(self.space)) * dx)

    def __call__(self):
        f = Function(self.space)
        b = f.vector().copy() # dolfin bug 889021
        for mat, x, alpha in self.matvecs:
            b_ = mat * self._as_vector(x)
            if alpha != 1:
                b_ *= alpha
            b += b_
        for vec in self.vecs:
            b += vec
        if self.forms:
            assemble(sum(self.forms), tensor=b, add_values=True, reset_sparsity=False)
        return b


