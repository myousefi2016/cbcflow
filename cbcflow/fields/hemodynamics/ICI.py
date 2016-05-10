#from cbcpost import Field, Threshold, Magnitude, Norm, DomainAvg, DomainSD, TimeAverage
#from dolfin import assemble, Constant, dx
from cbcpost import MetaField2, Field, SubFunction, ConstantField, Threshold, DomainAvg, TimeAverage, Magnitude
from cbcpost.utils.slice import create_slice
from dolfin import GenericFunction, Function, dot, project, Constant, assemble, dx
from numpy import dot as npdot

class Dot(MetaField2):
    def compute(self, get):
        u1 = get(self.valuename1)
        u2 = get(self.valuename2)

        if u1 == None or u2 == None:
            return

        if not (isinstance(u1, GenericFunction) or isinstance(u2, GenericFunction)):
            return npdot(u1,u2)

        if not isinstance(u1, GenericFunction):
            u1 = Constant(u1)
            u1,u2 = u2,u1
        if not isinstance(u2, GenericFunction):
            u2 = Constant(u2)

        if isinstance(u2, Function):
            u1,u2 = u2,u1

        assert isinstance(u1, Function)
        assert isinstance(u2, GenericFunction)

        if u1.value_rank() == u2.value_rank():
            if u1.value_rank() == 0:
                V = u1.function_space()
            else:
                V = u1.function_space().sub(0).collapse()
        elif u1.value_rank() > u2.value_rank():
            assert u2.value_rank() == 0
            V = u1.function_space()
            u1,u2 = u2,u1
        else:
            assert isinstance(u2, Function)
            assert u1.value_rank() == 0
            V = u2.function_space()

        N = max([u1.value_rank(), u2.value_rank()])

        if not hasattr(self, "u"):
            self.u = Function(V)

        if isinstance(u2, Function) and u1.function_space().dim() == u2.function_space().dim():
            self.u.vector()[:] = u1.vector().array()*u2.vector().array()
        elif u1.value_rank() == u2.value_rank():
            project(dot(u1,u2), function=self.u)
        else:
            assert u1.value_rank() == 0
            if isinstance(u1, Constant):
                self.u.vector()[:] = float(u1)*u2.vector().array()
            else:
                project(u1*u2, function=self.u)

        return self.u




class ICI(Field):
    def __init__(self, neck, pa_planes, *args, **kwargs):
        Field.__init__(self, *args, **kwargs)
        self.neck, self.necknormal = neck[0], neck[1]
        self.pa_planes = pa_planes

    @classmethod
    def default_params(cls):
        params = Field.default_params()
        params.update(use_timeaverage=False,
                      debug=False)
        return params

    def add_fields(self):
        Field.start_recording()

        params = self.params.copy_recursive()
        if not self.params.debug:
            params["save"] = False
            params["plot"] = False
        params.pop("debug")
        params.pop("use_timeaverage")
        params.pop("finalize")
        T0, T1 = self.params.start_time, self.params.end_time
        assert T0 != Field.default_params().start_time
        assert T1 != Field.default_params().end_time
        if self.params.use_timeaverage:
            u = TimeAverage("Velocity", params=params)
            velocity = u.name
        else:
            velocity = "Velocity"
        

        uneck = SubFunction(velocity, self.neck, params=params, label="neck")

        A = assemble(Constant(1)*dx(domain=self.neck))
        Qin = Dot(ConstantField(self.necknormal), uneck, params=params)

        #Qin.params.update(params)
        t = Threshold(Qin, ConstantField(0), dict(threshold_by="above"))
        t.params.update(params)
        t.name = "threshold_neck"

        Ain = A*DomainAvg(t)
        Ain.name = "Ain"
        Ain.params.update(params)
        Qin = A*DomainAvg(Dot(Qin,t, params=params))
        Qin.name = "Qin"
        Qin.params.update(params)
        
        Qpa = 0
        for i, (plane, n) in enumerate(self.pa_planes):
            upa = SubFunction(velocity, plane, params=params, label="pa_%d" %i)
            
            Ai = assemble(Constant(1)*dx(domain=plane))
            Q = Ai*DomainAvg(Dot(ConstantField(n), upa, params=params))
            Q.name = "Qpa%d" %i
            Q.params.update(params)
            Qpa += Magnitude(Q)
            Qpa.params.update(params)

        f = (Qin/Qpa)/(Ain/A)
        if not self.params.use_timeaverage:
            f = TimeAverage(f.name, params=params)
        
        self.valuename = f.name
        
        self.Qin = Qin.name
        self.Qpa = Qpa.name
        self.Ain = Ain.name
        self.A = A
        

        fields = Field.stop_recording()

        return fields

    def compute(self, get):
        Ain = get(self.Ain)
        A = self.A
        Qin = get(self.Qin)
        Qpa = get(self.Qpa)
        
        print "Ain: ", Ain
        print "A: ", A
        print "Qin: ", Qin
        print "Qpa: ", Qpa
        
        sci = get(self.valuename)
        #print "SCI: ", sci
        return sci