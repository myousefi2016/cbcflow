
from cbcflow.core.paramdict import ParamDict

def error(msg):
    print msg
    crash

class Parameterized(object):
    def __init__(self, params=None):
        cls = type(self) # Just to make it clear that merge_params is a class method
        self._params = cls.merge_params() if params is None else cls.merge_params(**params)

    @property
    def params(self):
        "Read-only access to instance parameters."
        return self._params

    @classmethod
    def class_params(cls):
        "Return two dicts (add, replace) with parameters to add for this class and parameters to override from subclass."
        raise NotImplementedError("A Parameterized subclass must implement class_params.")

    @classmethod
    def merge_params(cls, *args_override, **kwargs_override):
        "Merge parameters from subclasses and given overrides. Not for overriding in subclasses."
        classes = [c for c in reversed(cls.mro()) if hasattr(c, "class_params")]
        classes.remove(Parameterized)
        p = {}
        for c in classes:
            add, replace = c.class_params()

            for k in replace:
                if k not in p:
                    error("Trying to override non-existing parameter {0}.".format(k))
                p[k] = replace[k]
            #p.update(replace)

            for k in add:
                if k in p:
                    error("Trying to add an already existing parameter {0}.".format(k))
                p[k] = add[k]
            #p.update(add)

        for k in kwargs_override:
            if k not in p:
                error("Trying to override non-existing parameter {0}.".format(k))
            p[k] = kwargs_override[k]
        #p.update(kwargs_override)

        for override in args_override + (kwargs_override,):
            for k in override:
                if k not in p:
                    error("Trying to override non-existing parameter {0}.".format(k))
                p[k] = override[k]
            #p.update(override)

        return ParamDict(p)
