
from ..bases.PPField import PPField

class MetaPPField(PPField):
    def __init__(self, value, params=None, label=None):
        PPField.__init__(self, params, label)
        self.valuename = value.name if isinstance(value, PPField) else value

    @property
    def name(self):
        n = "%s_%s" % (self.__class__.__name__, self.valuename)
        if self.label: n += "_"+self.label
        return n

