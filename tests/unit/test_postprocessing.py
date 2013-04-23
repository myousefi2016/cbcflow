"""
Tests of the postprocessing framework in headflow.
"""

import unittest
from headflow import PostProcessorBase, PPFieldBase, ParamDict

class MockPostProcessor(PostProcessorBase):
    def __init__(self):
        PostProcessorBase.__init__(self)
        
    def print_data(self):
        for inst in self.list_all:
            print inst.get_data()
            
    def print_all_params(self):
        for inst in self.list_all:
            print inst.params


class MockA(PPFieldBase):
   
    def __init__(self, **kwargs):
        PPFieldBase.__init__(self, **kwargs)
    
    def update(self, u, p, t, timestep, problem):
        value = timestep+10
        self.set_data(t, timestep, value)
        
        

class MockB(PPFieldBase):
    def __init__(self, **kwargs):
        assert("parent" in kwargs.keys())
        PPFieldBase.__init__(self, **kwargs)
       
    def update(self, u, p, t, timestep, problem):
        # Get parent data
        parent_datadict = self.parent.get_data()
        parent_data = parent_datadict["data"]
        
        value = parent_data**0.5
        self.set_data(t, timestep, value)
        
    
class MockC(PPFieldBase):
    def __init__(self, **kwargs):
        assert("parent" in kwargs.keys())
        PPFieldBase.__init__(self, **kwargs)

    def update(self, u, p, t, timestep, problem):
        # Get parent data
        parent_datadict = self.parent.get_data()
        parent_data = parent_datadict["data"]
        
        value = parent_data**2
        self.set_data(t, timestep, value)
        
        
class TestPostProcessing(unittest.TestCase):
    def test_postprocessing(self):
        # Create postprocessor
        PP = MockPostProcessor()
                    
        # MockB and MockC needs a MockA instance as parent
        a = MockA(timeparams=ParamDict(start_timestep=5, end_timestep=15))
        b = MockB(parent=a, timeparams=ParamDict(start_timestep=1, end_timestep=10))
        c = MockC(parent=a, timeparams=ParamDict(start_timestep=5, end_timestep=18, step_frequency=4))
        
        # Add fields to postprocessor
        PP.add_field(a)
        PP.add_field(b)
        PP.add_field(c)
        
        from numpy import linspace
        t_range = linspace(0,1,21)
        
        # Dummy timeloop
        for timestep, t in enumerate(t_range):
            #print "######### Finished timestep %d (t=%f) ###########" %(timestep, t)
            PP.update_all(None, None, t, timestep, None)
            #PP.print_data()
        #PP.print_all_params()