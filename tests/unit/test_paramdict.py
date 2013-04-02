
import unittest

from headflow import ParamDict

class TestParamDict(unittest.TestCase):
    def test_init_from_dict(self):
        d = { 'a': 1, 'b': 3.14 }
        pd = ParamDict(d)
        self.assertEqual(len(pd), 2)
        self.assertEqual(pd.a, d['a'])
        self.assertEqual(pd.b, d['b'])
        self.assertEqual(pd['a'], d['a'])
        self.assertEqual(pd['b'], d['b'])

    def test_init_by_kwargs(self):
        pd = ParamDict(foo='hei', bar='argh')
        self.assertEqual(pd.foo, 'hei')
        self.assertEqual(pd.bar, 'argh')
        self.assertEqual(pd["foo"], 'hei')
        self.assertEqual(pd["bar"], 'argh')
        self.assertEqual(len(pd), 2)

    def test_init_by_sequence(self):
        keys = ('a', 'b', 'c')
        values = (1, 2, 3)
        items = tuple(zip(keys, values))
        pd = ParamDict(items)
        self.assertEqual(pd.a, 1)
        self.assertEqual(pd.b, 2)
        self.assertEqual(pd.c, 3)

    def test_add_params_after_init(self):
        pd = ParamDict()
        pd.a = 1
        pd.b = 2
        self.assertEqual(pd.a, 1)
        self.assertEqual(pd.b, 2)
        self.assertEqual(len(pd), 2)

    def test_shallow_iteration(self):
        keys = ('a', 'b', 'c')
        values = (1, 2, 3)
        items = tuple(zip(keys, values))
        pd = ParamDict(items)
        self.assertEqual(tuple(sorted(pd)), keys)
        self.assertEqual(tuple(sorted(pd.iterkeys())), keys)
        self.assertEqual(tuple(sorted(pd.keys())), keys)
        self.assertEqual(tuple(sorted(pd.iteritems())), items)
        self.assertEqual(tuple(sorted(pd.itervalues())), values)

    def create_multilevel_pd(self):
        pda1 = ParamDict(a=1)
        pdb1 = ParamDict(b=2)
        pdc1 = ParamDict(pa=pda1, pb=pdb1)
        pda2 = ParamDict(a=3)
        pdb2 = ParamDict(b=4)
        pdc2 = ParamDict(pa=pda2, pb=pdb2)
        pdd = ParamDict(pc1=pdc1, pc2=pdc2)
        return pdd
    
    def test_multilevel_access(self):
        pdd = self.create_multilevel_pd()
        self.assertEqual(pdd.pc1.pa.a, 1)
        self.assertEqual(pdd.pc1.pb.b, 2)
        self.assertEqual(pdd.pc2.pa.a, 3)
        self.assertEqual(pdd.pc2.pb.b, 4)

    def test_iterdeep_shallow_data(self):
        pd = ParamDict()
        deep = tuple(pd.iterdeep())
        self.assertEqual(deep, ())

        pd = ParamDict(a=3, b=4)
        deep = tuple(pd.iterdeep())
        self.assertEqual(deep, (('a',3), ('b',4)))

    def test_iterdeep_multilevel_data(self):
        pdd = self.create_multilevel_pd()
        deep = tuple(sorted(pdd.iterdeep()))
        items = ( ('pc1.pa.a', 1),
                  ('pc1.pb.b', 2),
                  ('pc2.pa.a', 3),
                  ('pc2.pb.b', 4), )
        self.assertEqual(deep, items)

    def test_shallow_copy(self):
        pd1 = ParamDict(a=3, b=4)
        pd2 = pd1.copy_recursive()
        pd1.a = 1
        pd2.a = 2
        pd2.b = 2
        pd1.b = 1
        self.assertEqual(pd1.a, 1)
        self.assertEqual(pd1.b, 1)
        self.assertEqual(pd2.a, 2)
        self.assertEqual(pd2.b, 2)

    def test_recursive_copy(self):
        pdcc1 = ParamDict(cca=30)
        pdcc2 = ParamDict(ccb=40)
        pdc1 = ParamDict(a=3, b=4, cc1=pdcc1, cc2=pdcc2)
        pdc2 = ParamDict(c=5, d=6)
        pd1 = ParamDict(c1=pdc1, c2=pdc2)
        pd2 = pd1.copy_recursive()

        self.assertEqual(pd1.c2.d, 6)
        self.assertEqual(pd2.c2.d, 6)
        pd1.c2.d = 7
        self.assertEqual(pd1.c2.d, 7)
        self.assertEqual(pd2.c2.d, 6)
        pd2.c2.d = 8
        self.assertEqual(pd1.c2.d, 7)
        self.assertEqual(pd2.c2.d, 8)

        self.assertEqual(pd1.c1.cc2.ccb, 40)
        self.assertEqual(pd2.c1.cc2.ccb, 40)
        pd2.c1.cc2.ccb = 50
        self.assertEqual(pd1.c1.cc2.ccb, 40)
        self.assertEqual(pd2.c1.cc2.ccb, 50)
        pd1.c1.cc2.ccb = 60
        self.assertEqual(pd1.c1.cc2.ccb, 60)
        self.assertEqual(pd2.c1.cc2.ccb, 50)

    def test_shallow_update(self):
        pd1 = ParamDict(a=3, b=4)
        pd2 = ParamDict(b=14, c=15)
        pd1orig = pd1.copy()
        pd1.update_shallow(pd2)
        self.assertTrue(all(k in pd1 for k in pd1orig))
        self.assertTrue(all(k in pd1 for k in pd2))
        self.assertTrue(all(pd1[k] == pd2[k] for k in pd2))
        self.assertTrue(all(pd1[k] == pd1orig[k] for k in pd1orig if not k in pd2))

    def test_recursive_update(self):
        # Build multilevel test data
        pdcc1 = ParamDict(cca=30)
        pdcc2 = ParamDict(ccb=40)
        pdc1 = ParamDict(a=3, b=4, cc1=pdcc1, cc2=pdcc2)
        pdc2 = ParamDict(c=5, d=6)
        pdorig = ParamDict(c1=pdc1, c2=pdc2, v=7)

        # Build alternative multilevel test data
        apdcc1 = ParamDict(cca=31)
        apdcc2 = ParamDict(ccb=41)
        apdc1 = ParamDict(a=5, b=8, cc1=apdcc1, cc2=apdcc2)
        apdc2 = ParamDict(c=7, d=9)
        apdorig = ParamDict(c1=apdc1, c2=apdc2, v=2)

        self.assertNotEqual(pdorig, apdorig)

        # Supply a single item
        pd = pdorig.copy_recursive()
        self.assertEqual(pd.v, 7)
        pd.update_recursive(v=9)
        self.assertEqual(pd.v, 9)
        pd.update_recursive({'v':10})
        self.assertEqual(pd.v, 10)
        pd.update_recursive(ParamDict(v=11))
        self.assertEqual(pd.v, 11)

        # Supply multiple items for child
        pd = pdorig.copy_recursive()
        self.assertEqual(pd.c1.a, 3)
        self.assertEqual(pd.c1.b, 4)
        if 0:
            pd.update_recursive(c1={'a':11, 'b':22})
            self.assertEqual(pd.c1.a, 11)
            self.assertEqual(pd.c1.b, 22)
        pd.update_recursive(c1=ParamDict(a=13, b=25))
        self.assertEqual(pd.c1.a, 13)
        self.assertEqual(pd.c1.b, 25)

        # Supply a full multilevel paramdict
        pd = pdorig.copy_recursive()
        self.assertEqual(pd, pdorig)
        self.assertNotEqual(pd, apdorig)
        pd.update_recursive(apdorig)
        self.assertNotEqual(pd, pdorig)
        self.assertEqual(pd, apdorig)

    def test_pickle_protocol(self):
        pass # TODO

    def test_repr_rendering(self):
        pass # TODO

    def test_str_rendering(self):
        pdid = ParamDict(a=1, b=2)
        pdin = ParamDict(a=1, b=3)
        pdout = ParamDict(a=1, b=4)
        record = ParamDict(identity=pdid,
                           input=pdin,
                           output=pdout)
        s = str(record)
        self.assertEqual(record["identity"], pdid)

    def test_arg_rendering(self):
        pass # TODO

    def test_arg_parsing(self):
        pass # TODO
