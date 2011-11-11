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
            self.matvecs.append((A, x, 1))
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
                self.matvecs.append((A, x, -1))
                return self
        raise RuntimeError, "Try '+=' instead"

    def _as_vector(self, x):
        if isinstance(x, GenericVector):
            return x
        if isinstance(x, Function):
            return x.vector()
        f = interpolate(x, self.space)
        v = f.vector()
        v._dummy = f # dolfin bug 889021
        return v

    def __call__(self):
        f = Function(self.space)
        b = f.vector()
        b._dummy = f # dolfin bug 889021
        for mat, x, alpha in self.matvecs:
            b_ = mat * self._as_vector(x)
            if alpha != 1:
                b_ *= alpha
            b += b_
        for vec in self.vecs:
            b += vec
        for form in self.forms:
            assemble(form, tensor=b, add_values=True, reset_sparsity=False)
        return b

