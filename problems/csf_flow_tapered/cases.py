from EllipticPipeMesh import EllipticPipeMesh

class casebase(object):
    def __init__(self):
        self.case = EllipticPipeMesh(self.mesh_name,
                                     self.a1, self.b1,
                                     self.a2, self.b2,
                                     self.z_min, self.z_max,
                                     self.x_index, self.y_index, self.z_index,
                                     a1_i=self.a1_i, b1_i=self.b1_i,
                                     a2_i=self.a2_i, b2_i=self.b2_i)

    def get_sub_domains(self):
        return self.case.get_sub_domains()

    def get_mesh(self):
        return self.case.get_mesh()

    def get_z_index(self):
        return self.z_index

    def get_contour(self):
        return self.case.get_contour()

    def get_top(self):
        return self.case.get_top()

    def get_bottom(self):
        return self.case.get_bottom()

    # default values, may be overridden in subclass
    a1 = 0.75
    a2 = 0.75
    b1 = 0.75
    b2 = 0.75
    z_min = -5.0
    z_max = 5.0
    x_index = 0
    y_index = 1
    z_index = 2
    a1_i = 0.5
    a2_i = 0.5
    b1_i = 0.5
    b2_i = 0.5

class caseStenosis40(casebase):
    mesh_name = "data/meshes/chiari/csf_block40.xml"

class stenosis40_refined(casebase):
    mesh_name = "data/meshes/chiari/block40_ref.xml"

class caseStenosis60(casebase):
    mesh_name = "data/meshes/chiari/csf_block60.xml"

class stenosis60_refined(casebase):
    mesh_name = "data/meshes/chiari/block60_ref.xml"

class stenosis75_refined(casebase):
    mesh_name = "data/meshes/chiari/block75_ref3.xml"

class stenosis85_ref(casebase):
    #mesh_name = "data/meshes/chiari/stenosis85_ref.xml.gz"
    mesh_name = "data/meshes/chiari/stenosis85_ng_ref2.xml.gz"

class straight_nerves(casebase):
    mesh_name = "data/meshes/chiari/straight_nerves_bl.xml.gz"
    z_min = -4.9718099999999996
    z_max = 4.9809900000000003
