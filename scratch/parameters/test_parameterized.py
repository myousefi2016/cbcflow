
from parameterized import Parameterized

def test_parameterized_classes():
    class Base(Parameterized):
        def __init__(self, params=None):
            Parameterized.__init__(self, params)

        @classmethod
        def class_params(cls):
            add = {"base": 0}
            replace = {}
            return add, replace

    class Sub(Base):
        def __init__(self, params=None):
            Base.__init__(self, params)

        @classmethod
        def class_params(cls):
            add = {"sub": 11}
            replace = {"base": 10}
            return add, replace

    class Sub2(Sub):
        def __init__(self, params=None):
            Sub.__init__(self, params)

        @classmethod
        def class_params(cls):
            add = {"sub2": 22}
            replace = {"base": 20, "sub":21}
            return add, replace

    assert Base.merge_params() == {"base": 0}
    assert Base.merge_params(base=2) == {"base": 2}
    assert Base({"base":2}).params == {"base": 2}
    assert Base().params == {"base": 0}
    assert Base().params.base == 0

    assert Sub.merge_params() == {"base": 10, "sub": 11}
    assert Sub.merge_params(sub=4) == {"base": 10, "sub": 4}
    assert Sub.merge_params(base=3, sub=4) == {"base": 3, "sub": 4}
    assert Sub().params == {"base": 10, "sub": 11}

    assert Sub2.merge_params() == {"base": 20, "sub": 21, "sub2": 22}
    assert Sub2.merge_params(sub2=4) == {"base": 20, "sub": 21, "sub2": 4}
    assert Sub2.merge_params(base=3, sub=5, sub2=4) == {"base": 3, "sub": 5, "sub2": 4}
    assert Sub2().params == {"base": 20, "sub": 21, "sub2": 22}
    assert Sub2().params.sub == 21
