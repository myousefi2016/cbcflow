import sys
sys.path.insert(0,"../../site-packages")

from headflow import Pouseille, Womersley
from headflow.dol import Function, VectorFunctionSpace, Mesh, Expression, DirichletBC

mesh = Mesh("cylinder_4k.xml.gz")

x = [0.0, 0.2, 0.4, 0.6, 0.8]
y = [1, 5, 3, 2, 1]

coeffs = zip(x,y)

V = VectorFunctionSpace(mesh, "CG", 1)
u = Function(V)

def test_bcs(bcs):
    # TODO: Expand
    for i, bc in enumerate(bcs):
        assert(isinstance(bc, Expression))
        dbc = DirichletBC(V.sub(i), bc, 1)
        dbc.apply(u.vector())
        
test_bcs(Pouseille(coeffs, mesh, 1))
test_bcs(Womersley(coeffs, mesh, 1, 4.0))


