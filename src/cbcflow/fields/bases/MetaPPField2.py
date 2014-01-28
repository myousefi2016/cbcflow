
from ..bases.PPField import PPField

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
