
from dolfin import *

def as_point(p):
    if isinstance(p, Point):
        return p
    else:
        return Point(*p)

def as_point_value(p):
    if isinstance(p, Point):
        return (p.x(), p.y())
    else:
        return p

def as_points(point_values):
    return [as_point(p) for p in point_values]

def as_point_values(points):
    return [as_point_value(p) for p in points]

def translate_points(offset, points):
    print offset
    offset = as_point(offset)
    return [as_point(p) + offset for p in points]

def scale_point(factor, point):
    point = as_point_value(point)
    return as_point((factor*point[0], factor*point[1]))

def scale_points(factor, points):
    return [scale_point(factor, p) for p in points]

def rotate_points(theta, points):
    axis = Point(0.0, 0.0, 1.0)
    return [as_point(p).rotate(axis, theta) for p in points]

def rectangle_points(ap, bp):
    a0, a1 = as_point_values(ap)
    b0, b1 = as_point_values(bp)
    return as_points([(a0, a1), (b0, a1), (b0, b1), (a0, b1)])

def unit_square_points():
    return rectangle_points((0.0, 0.0), (1.0, 1.0))


def make_mesh_from_points(points, res=16):
    polygon = Polygon(as_points(points))
    mesh = Mesh(polygon, res)
    return mesh, None

def make_diamond1(res=16):
    mesh_res = 2*res

    points = as_points([(0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0)])

    polygon = Polygon(points)
    mesh = Mesh(polygon, res)
    return mesh, None

def make_diamond(res=16):

    theta = pi/4
    factor = 0.5**0.5
    offset = (0.0, 1.0)

    points = rotate_points(theta, points)
    points = scale_points(factor, points)
    points = translate_points(offset, points)

    polygon = Polygon(points)

    mesh = Mesh(polygon, res)
    return mesh, None

def make_angled_channel(res=16):
    mesh_res = 2*res

    origin = (0.0, 0.0)
    angle = pi/4
    radius = 0.5
    length = 10.0

    points = rectangle_points((0.0, -radius), (length, +radius))
    points = rotate_points(angle, points)
    points = translate_points(origin, points)

    polygon = Polygon(points)

    mesh = Mesh(polygon, mesh_res)
    return mesh, None

def make_channel(res=16):
    mesh_res = 2*res

    channel = Rectangle(0.0, -0.5, 10.0, +0.5)
    geometry = channel

    mesh = Mesh(geometry, mesh_res)
    return mesh, None

def make_channel_with_aneurysm(res=16):
    cres = res
    mesh_res = 2*res

    channel = Rectangle(0.0, -0.5, 10.0, +0.5)
    aneurysm = Circle(5.0, 0.6, 0.7, cres)
    geometry = channel + aneurysm

    mesh = Mesh(geometry, mesh_res)
    return mesh, None

def make_t_section(res=16):
    cres = res
    mesh_res = 2*res

    inlet = Rectangle(-0.5, 0.0, +0.5, 10.0)
    outlets = Rectangle(-5.0, 9.0, +5.0, 10.0)
    geometry = inlet + outlets

    mesh = Mesh(geometry, mesh_res)
    return mesh, None

def make_t_section_with_aneurysm(res=16):
    cres = res
    mesh_res = 2*res

    inlet = Rectangle(-0.5, 0.0, +0.5, 10.0)
    outlets = Rectangle(-5.0, 9.0, +5.0, 10.0)
    aneurysm = Circle(0.0, 10.1, 0.7, cres)
    geometry = inlet + outlets + aneurysm

    mesh = Mesh(geometry, mesh_res)
    return mesh, None

def make_bifurcation_with_aneurysm(res=16):
    # TODO: Make this a reusable function with parameterization as arguments

    # Mesh and circle resolutions
    cres = res
    pres = res
    mesh_res = 2*res

    # Parameterization of geometry
    midx, midy = 0.0, 0.0
    radius = 0.5
    inlen = 8.0

    left_outlen = 6.0
    left_angle = pi/4
    left_relrad = 0.5**0.5 * 1.01

    right_outlen = 6.0
    right_angle = pi/4
    right_relrad = 0.5**0.5 * 1.01

    arad = 1.5
    aoff = 2.5
    aside = 0.2

    # (origin, angle, length, relradius) for each channel segment
    segments = [
        ((midx, midy), -pi/2, inlen, radius),
        ((midx, midy), pi/2 + left_angle, left_outlen, left_relrad*radius),
        ((midx, midy), pi/2 - right_angle, right_outlen, right_relrad*radius),
        ]

    # Assert that we don't get gaps in the angles between the inlet and outlet channels
    assert left_relrad / cos(pi/2-left_angle) > 1.0, "Gap between inlet and left outlet channel."
    assert right_relrad / cos(pi/2-right_angle) > 1.0, "Gap between inlet and right outlet channel."

    # Create channel segment polygons
    channels = []
    for origin, angle, length, r in segments:
        points = rectangle_points((0.0, -r), (length, +r))
        points = rotate_points(angle, points)
        points = translate_points(origin, points)
        channels += [Polygon(points)]

    # Create aneurysm
    ax, ay = aside*radius, aoff*radius
    ar = arad*radius
    aneurysm = Circle(ax, ay, ar, cres)

    # Create patches to avoid holes
    patch1 = Circle(midx, midy, radius, pres)
    patch2 = Circle(ax, ay-ar, 0.5*radius, pres)

    # Merge geometries
    geometry = channels[0]
    for channel in channels[1:]:
        geometry = geometry + channel
    geometry = geometry + aneurysm
    geometry = geometry + patch1
    geometry = geometry + patch2

    # Discretize
    mesh = Mesh(geometry, mesh_res)

    # Mark boundaries
    class RadiusDomain(SubDomain):
        def __init__(self, x, r):
            SubDomain.__init__(self)
            self.x = x
            self.r2 = r**2 * (1.0+10*DOLFIN_EPS)
        def inside(self, x, on_boundary):
            r2 = (x[0]-self.x[0])**2 + (x[1]-self.x[1])**2
            return on_boundary and r2 <= self.r2
    facet_domains = FacetFunction("size_t", mesh)
    facet_domains.set_all(2+len(segments))
    DomainBoundary().mark(facet_domains, 0)
    RadiusDomain((ax, ay), ar).mark(facet_domains, 1)
    for k, seg in enumerate(segments):
        origin, angle, length, r = seg
        x = origin[0] + length*cos(angle)
        y = origin[1] + length*sin(angle)
        RadiusDomain((x,y), r).mark(facet_domains, 2+k)

    return mesh, facet_domains

#make_mesh = make_channel
#make_mesh = make_channel_with_aneurysm
#make_mesh = make_t_section
#make_mesh = make_t_section_with_aneurysm
#make_mesh = make_bifurcation
#make_mesh = make_diamond
#make_mesh = make_angled_channel
make_mesh = make_bifurcation_with_aneurysm

def create_meshes():
    resolutions = (8, 16, 32, 64)
    meshes = []
    for k, res in enumerate(resolutions):
        mesh, facet_domains = make_mesh(res)
        meshes.append((res, mesh, facet_domains))
    return meshes

def store_meshes(meshes):
    for res, mesh, facet_domains in meshes:
        File("mesh_r%d.xml.gz" % res) << mesh
        if facet_domains is not None:
            File("facet_domains_r%d.xml.gz" % res) << facet_domains

def count_meshdims(meshes):
    for res, mesh, facet_domains in meshes:
        U = FunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)
        R = FunctionSpace(mesh, "CR", 1)
        P = FunctionSpace(mesh, "DG", 0)
        print
        print "Sizes for mesh with res = %d:" % (res,)
        print "c,v:", mesh.num_cells(), mesh.num_vertices()
        print "u,q,uq:", U.dim(), Q.dim(), U.dim()*3 + Q.dim()
        print "r,p,rp:", R.dim(), P.dim(), R.dim()*3 + P.dim()
        print "Assembled volume of mesh:"
        print assemble(Constant(1.0)*dx(), mesh=mesh)
        if facet_domains is not None:
            print "Assembled areas of boundaries:"
            for k in range(6):
                print k, assemble(Constant(1.0)*ds[facet_domains](k), mesh=mesh)

def plot_meshes(meshes):
    for res, mesh, facet_domains in meshes:
        plot(mesh, title="mesh %d" % res)
        if facet_domains is not None:
            plot(facet_domains, title="facet domains %d" % res)
    interactive()

if __name__ == "__main__":
    meshes = create_meshes()
    count_meshdims(meshes)
    store_meshes(meshes)
    plot_meshes(meshes)
