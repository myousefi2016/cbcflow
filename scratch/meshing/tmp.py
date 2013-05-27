from dolfin import *

def rotate_points(points, theta):
    return [Point(p).rotate(Point(0.0,0.0,1.0), theta) for p in points]

def rectangle_points(x0, x1, y0, y1):
    return [(x0, x1), (y0, x1), (y0, y1), (x0, y1)]

def unit_square_points():
    return rectangle_points(0.0, 0.0, 1.0, 1.0)

points = unit_square_points()

points = [Point(*p) for p in points]
#points = [(p.x(), p.y(), p.z()) for p in points]

res = 16
polygon = Polygon(points)
mesh = Mesh(polygon, res)
plot(mesh)
interactive()

