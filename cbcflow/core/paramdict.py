from __future__ import division
__author__ = "Martin Alnaes <martinal@simula.no>"
__date__ = "2013-04-26"
__copyright__ = "Copyright (C) 2011-2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

import re

class ParamDict(dict):
    def __init__(self, *args, **kwargs):
        dict.__init__(self)
        self._keys = []
        if args:
            arg, = args
            for item in arg:
                if isinstance(item, tuple):
                    k, v = item
                else:
                    k, v = item, arg[item]
                self[k] = v
        if kwargs:
            self.update_recursive(kwargs)
            self._keys = sorted(set(self._keys) | set(kwargs))

    # --- Recursive ParamDict aware copy and update functions

    def copy_recursive(self):
        "Copy ParamDict hierarchy recursively, using copy.deepcopy() to copy values."
        import copy
        keys = list(self.iterkeys())
        items = []
        for k in keys:
            v = self[k]
            if isinstance(v, ParamDict):
                v2 = v.copy_recursive()
            else:
                v2 = copy.deepcopy(v)
            items.append((k,v2))
        return ParamDict(items)

    def replace_shallow(self, params=None, **kwparams):
        "Perform a shallow update where no new keys are allowed."
        if params:
            unknown = set(params.iterkeys()) - set(self.iterkeys())
            if unknown:
                raise RuntimeError("Trying to replace non-existing entries: %s" % (sorted(unknown),))
            for k, v in params.iteritems():
                self[k] = v
        if kwparams:
            self.replace_shallow(kwparams)

    def replace_recursive(self, params=None, **kwparams):
        "Perform a recursive update where no new keys are allowed."
        def handle(k, v):
            if k not in self:
                raise RuntimeError("Trying to replace non-existing entry: %s" % (k,))
            if isinstance(v, ParamDict):
                # If it's a ParamDict, recurse
                self[k].replace_recursive(v)
            else:
                # Otherwise abort recursion
                self[k] = v
        if params:
            for k, v in params.iteritems():
                handle(k, v)
        for k, v in kwparams.iteritems():
            handle(k, v)

    def update_shallow(self, params=None, **kwparams):
        "Perform a shallow update, allowing new keys to be introduced."
        if params:
            for k, v in params.iteritems():
                self[k] = v
        if kwparams:
            self.update_shallow(kwparams)

    def update_recursive(self, params=None, **kwparams):
        "Perform a recursive update, allowing new keys to be introduced."
        def handle(k, v):
            if isinstance(v, ParamDict):
                # If it's a ParamDict, recurse
                pd = ParamDict()
                pd.update_recursive(v)
                self[k] = pd
            else:
                # Otherwise abort recursion
                self[k] = v
        if params:
            for k, v in params.iteritems():
                handle(k, v)
        for k, v in kwparams.iteritems():
            handle(k, v)

    update = update_recursive
    replace = replace_recursive

    # --- Attribute access

    def __getitem__(self, name):
        return dict.__getitem__(self, name)

    def __setitem__(self, name, value):
        "Insert item with dict notation, allowing a new key to be added."
        if name not in self:
            assert(isinstance(name, str))
            self._keys.append(name)
        return dict.__setitem__(self, name, value)

    def __delitem__(self, name): # TODO: Add tests for this
        if name in self:
            assert(isinstance(name, str))
            self._keys.remove(name)
        return dict.__delitem__(self, name)

    def __getattr__(self, name):
        if name.startswith("_"):
            return self.__dict__[name]
        else:
            return self[name]

    def __setattr__(self, name, value):
        "Insert item with attribute notation, only allows changing a value with existing key."
        if name.startswith("_"):
            self.__dict__[name] = value
        else:
            if name not in self:
                raise RuntimeError("Trying to update non-existing entry: %s" % (name,))
            self[name] = value

    def pop(self, name, default=None):
        ''' Returns Paramdict[name] if the key exists. If the key does not exist the default value is returned. '''
        if self.has_key(name):
            v = self[name]
            del self[name]
            return v
        else:
            return default

    # --- Pickling and shelving

    def __getstate__(self):
        return (self._keys, list(dict.iteritems(self)))

    def __setstate__(self, state):
        (self._keys, data) = state
        assert len(self._keys) == len(data)
        dict.clear(self)
        dict.update(self, data)

    # --- String rendering

    def __repr__(self):
        return "ParamDict([%s])" % ", ".join("(%r, %r)" % (k, self[k]) for k in self._keys)

    def __str__(self):
        return '\n'.join(self._str())

    def _str(self, level=0):
        indent = (level+1)*4*" "
        if level == 0:
            lines = ["Parameters:"]
        else:
            lines = []

        for k in self._keys:
            v = self[k]
            if isinstance(v, ParamDict):
                lines.append("%s%s =" % (indent, k))
                lines.extend(v._str(level+1))
            else:
                lines.append("%s%s = %r" % (indent, k, v))

        return lines

    # TODO: Add .json format support!

    # --- Iteration

    def __iter__(self):
        return iter(self._keys)

    def iteritems(self):
        return ((k, self[k]) for k in self._keys)

    def items(self):
        return list(self.iteritems())

    def iterkeys(self):
        return iter(self._keys)

    def keys(self):
        return list(self._keys)

    def iterdeep(self):
        "Iterate recursively over all parameter items."
        for k, v in self.iteritems():
            if isinstance(v, ParamDict):
                for sk, sv in v.iterdeep():
                    yield ("%s.%s" % (k, sk), sv)
            else:
                yield (k, v)

    # --- Commandline argument translation (should perhaps be placed outside class?)

    def arg_assign(self, name, value):
        subs = name.split('.')
        subs, vname = subs[:-1], subs[-1]
        p = self
        for s in subs:
            p = p[s]
        p[vname] = eval(value)

    def parse_args(self, args):
        m = re.findall(r'([^ =]+)=([^ ]+)', args)
        for k, v in m:
            self.arg_assign(k, v)

    def render_args(self):
        return "  ".join("%s=%r" % (k,v) for k, v in self.iterdeep())