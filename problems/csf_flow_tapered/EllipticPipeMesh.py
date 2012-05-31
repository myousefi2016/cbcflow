from dolfin import SubDomain, Mesh, MeshFunction
from problems.problembase import ProblemBase

class MeshInfo(object):

	def __init__(self, mesh_name):
		self.mesh = Mesh(ProblemBase(None).retrieve(mesh_name))
		self.mesh.order()
		self.mesh.init(2)
        # Load sub domain markers
		self.sub_domains =  MeshFunction("uint", self.mesh, self.mesh.topology().dim() - 1)

	def get_mesh(self):
		return self.mesh

	def get_sub_domains(self):
		return self.sub_domains

# class to define mesh boundaries of flow in elliptic pipes 

# a1, b1, a2, b2 describe the vectors in the directons of the axis indicated by
# x_index (for a) and y_index (for b). The number indicates the inflow (1) and
# the outflow(2); z_index is the index of the axis pointing along the tube,
# z_max and z_min define the tubes starting and ending point in tubular
# direction; Moved center points are indicated by the direction (x/y), the
# position (1:inflow: 2: outflow)

# for flow between en inner and an outer tube, ...

# Susanne Hentschel 13 May 2009


def on_ellipse(x, a, b, x_index, y_index, x_move=0, y_move=0):
	# move the coordinates into so that (0,0) is the center of the ellipse
	x1 = x[x_index] - x_move
	x2 = x[y_index] - y_move
	return bool( abs((x1/a)**2 + (x2/b)**2 - 1.0 ) < 2*10**(-2) )

class Top(SubDomain):	#bc for top
	def __init__(self, a2_o, a2_i, b2_o, b2_i,  x_index, y_index, z_index, z_max, x2_o_move=0, y2_o_move=0, x2_i_move=0, y2_i_move=0):
		SubDomain.__init__(self)
		self.x_index = x_index
		self.y_index = y_index
		self.a2_o = a2_o
		self.a2_i = a2_i
		self.b2_o = b2_o
		self.b2_i = b2_i
		self.z_index = z_index
		self.z_max = z_max
		self.x2_o_move = x2_o_move
		self.x2_i_move = x2_i_move
		self.y2_o_move = y2_o_move
		self.y2_i_move = y2_i_move
	def inside(self, x, on_boundary):
		return bool(on_boundary and abs(x[self.z_index] - self.z_max) < 10**(-3)  and not on_ellipse(x, self.a2_o, self.b2_o, self.x_index, self.y_index, self.x2_o_move, self.y2_o_move ) and not on_ellipse(x, self.a2_i, self.b2_i, self.x_index, self.y_index, self.x2_i_move, self.y2_i_move ) )

class Bottom(SubDomain):	# bc for bottom
	def __init__(self, a1_o, a1_i, b1_o, b1_i, x_index, y_index, z_index, z_min, x1_o_move=0, y1_o_move=0, x1_i_move=0, y1_i_move=0):
		SubDomain.__init__(self)
		self.x_index = x_index
		self.y_index = y_index
		self.z_index = z_index
		self.z_min = z_min
		self.a1_o = a1_o
		self.a1_i = a1_i
		self.b1_o = b1_o
		self.b1_i = b1_i
		self.x1_o_move = x1_o_move
		self.x1_i_move = x1_i_move
		self.y1_o_move = y1_o_move
		self.y1_i_move = y1_i_move
	def inside(self, x, on_boundary):
		return bool(on_boundary and abs(x[self.z_index] - self.z_min)< 10**(-3) and not on_ellipse(x, self.a1_o, self.b1_o, self.x_index, self.y_index, self.x1_o_move, self.y1_o_move ) and not on_ellipse(x, self.a1_i, self.b1_i, self.x_index, self.y_index, self.x1_i_move, self.y1_i_move ))

class Contour(SubDomain):	# bc for rest
	def __init__(self, x_index, y_index, top, bottom):
		SubDomain.__init__(self)
		self.x_index = x_index
		self.y_index = y_index
		self.top = top
		self.bottom = bottom
	def inside(self, x, on_boundary):
		return bool(on_boundary and not self.top.inside(x, on_boundary) and not self.bottom.inside(x, on_boundary) )





class EllipticPipeMesh(MeshInfo):
	def __init__(self, mesh_name, a1, b1, a2, b2, z_min, z_max, x_index, y_index, z_index, x1_move=0, y1_move=0, x2_move=0, y2_move=0, a1_i=0, b1_i=0, a2_i=0, b2_i=0, x1_i_move=0, y1_i_move=0, x2_i_move=0, y2_i_move=0): 
		MeshInfo.__init__(self, mesh_name)
	#	self.boundary_name = boundary_name
		self.a1 = a1
		self.a2 = a2
		self.b1 = b1
		self.b2 = b2
		self.z_min = z_min
		self.z_max = z_max
		self.x_index = x_index
		self.y_index = y_index
		self.z_index = z_index
		self.x1_move = x1_move
		self.x2_move = x2_move
		self.y1_move = y1_move
		self.y2_move = y2_move
		self.a1_i = a1_i
		self.a2_i = a2_i
		self.b1_i = b1_i
		self.b2_i = b2_i
		self.x1_i_move = x1_i_move
		self.x2_i_move = x2_i_move
		self.y1_i_move = y1_i_move
		self.y2_i_move = y2_i_move


		
		# Mark all facets as sub domain 3
		for i in range(self.sub_domains.size()):
#			self.sub_domains.set(i, 3)
			self.sub_domains[i]=3
			
		
		self.top = Top(self.a2, self.a2_i, self.b2, self.b2_i, self.x_index, self.y_index, self.z_index, self.z_max, x2_o_move=self.x2_move, y2_o_move=self.y2_move, x2_i_move=self.x2_i_move, y2_i_move=self.y2_i_move )
		self.bottom = Bottom(self.a1, self.a1_i, self.b1, self.b1_i, self.x_index, self.y_index, self.z_index, self.z_min, x1_o_move=self.x1_move, y1_o_move=self.y1_move, x1_i_move=self.x1_i_move, y1_i_move=self.y1_i_move  )
		self.contour = Contour(self.x_index, self.y_index, self.top, self.bottom)




		self.contour.mark(self.sub_domains, 0)
		self.top.mark(self.sub_domains, 1)
		self.bottom.mark(self.sub_domains, 2)



	def get_contour(self):
		return self.contour

	def get_top(self):
		return self.top

	def get_bottom(self):
		return self.bottom

"""

	def get_contour(self):
		return self.contour

	def get_top(self):
		return self.top

	def get_bottom(self):
		return self.bottom


"""

