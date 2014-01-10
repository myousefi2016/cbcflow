
from .PPField import PPField

class MetaPPField(PPField):
    def __init__(self, value, params=None, label=None):
        PPField.__init__(self, params, label)
        self.valuename = value.name if isinstance(value, PPField) else value

    @property
    def name(self):
        n = "%s_%s" % (self.__class__.__name__, self.valuename)
        if self.label: n += "_"+self.label
        return n

class MetaPPField2(PPField):
    def __init__(self, value1, value2, params=None, label=None):
        PPField.__init__(self, params, label)
        self.valuename1 = value1.name if isinstance(value1, PPField) else value1
        self.valuename2 = value2.name if isinstance(value2, PPField) else value2

    @property
    def name(self):
        n = "%s_%s_%s" % (self.__class__.__name__, self.valuename1, self.valuename2)
        if self.label: n += "_"+self.label
            
        return "%s_%s_%s" % (self.__class__.__name__, self.valuename1, self.valuename2)
