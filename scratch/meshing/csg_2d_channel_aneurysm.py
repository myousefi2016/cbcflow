
from dolfin import *

def make_channel(res=16):
    mesh_res = 2*res

    channel = Rectangle(0.0, -0.5, 10.0, +0.5)
    geometry = channel

    mesh = Mesh(geometry, mesh_res)
    return mesh

def make_channel_with_aneurysm(res=16):
    cres = res
    mesh_res = 2*res

    channel = Rectangle(0.0, -0.5, 10.0, +0.5)
    aneurysm = Circle(5.0, 0.6, 0.7, cres)
    geometry = channel + aneurysm

    mesh = Mesh(geometry, mesh_res)
    return mesh

def make_t_section(res=16):
    cres = res
    mesh_res = 2*res

    inlet = Rectangle(-0.5, 0.0, +0.5, 10.0)
    outlets = Rectangle(-5.0, 9.0, +5.0, 10.0)
    geometry = inlet + outlets

    mesh = Mesh(geometry, mesh_res)
    return mesh

def make_t_section_with_aneurysm(res=16):
    cres = res
    mesh_res = 2*res

    inlet = Rectangle(-0.5, 0.0, +0.5, 10.0)
    outlets = Rectangle(-5.0, 9.0, +5.0, 10.0)
    aneurysm = Circle(0.0, 10.1, 0.7, cres)
    geometry = inlet + outlets + aneurysm

    mesh = Mesh(geometry, mesh_res)
    return mesh

def make_bifurcation(res=16):
    mesh_res = 2*res

    points = [
        (0.0,  0.0),
        (10.0, 0.0),
        (10.0, 2.0),
        (8.0,  2.0),
        (7.5,  1.0),
        (2.5,  1.0),
        (2.0,  4.0),
        (0.0,  4.0),
        (0.0,  0.0),
        ]
    points = [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0)]
    points = reversed(points)
    points = [Point(*p) for p in points]
    p = Polygon(points)
    geometry = p

    mesh = Mesh(geometry, mesh_res)
    return mesh

#make_mesh = make_channel
#make_mesh = make_channel_with_aneurysm
#make_mesh = make_t_section
#make_mesh = make_t_section_with_aneurysm
make_mesh = make_bifurcation

def create_meshes():
    resolutions = (8, 16, 32, 64)
    resolutions = (32,)
    meshes = []
    for k, n in enumerate(resolutions):
        mesh = make_mesh(n)
        meshes.append(mesh)

        U = FunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)
        R = FunctionSpace(mesh, "CR", 1)
        P = FunctionSpace(mesh, "DG", 0)
        print
        print "Sizes for mesh %d, resolution %d:" % (k, n)
        print mesh.num_cells(), mesh.num_vertices()
        print U.dim(), Q.dim(), U.dim()*3 + Q.dim()
        print R.dim(), P.dim(), R.dim()*3 + P.dim()
    return meshes

def plot_meshes(meshes):
    for k, mesh in enumerate(meshes):
        plot(mesh, title="mesh %d" % k)
    interactive()

if __name__ == "__main__":
    plot_meshes(create_meshes())

