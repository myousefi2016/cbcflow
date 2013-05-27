
from dolfin import *

def make_pipe_with_aneurysm(cres=16):
    mres = 2*cres

    r = Rectangle(0.0, -0.5, 10.0, +0.5)
    c = Circle(5.0, 0.6, 0.7, cres)
    g = r + c
    mesh = Mesh(g, mres)
    return mesh

meshes = []
for n in (8, 16, 32, 64):
    mesh = make_pipe_with_aneurysm(n)
    meshes.append(mesh)
    U = FunctionSpace(mesh, "CG", 2)
    Q = FunctionSpace(mesh, "CG", 1)
    R = FunctionSpace(mesh, "CR", 1)
    P = FunctionSpace(mesh, "DG", 0)
    print
    print "Sizes for mesh resolution %d:" % n
    print mesh.num_cells(), mesh.num_vertices()
    print U.dim(), Q.dim(), U.dim()*3 + Q.dim()
    print R.dim(), P.dim(), R.dim()*3 + P.dim()
    plot(mesh, title="mesh %d" % n)

interactive()
