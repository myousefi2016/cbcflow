import sys
sys.path.insert(0,"../../site-packages")

from dolfin import *
set_log_level(100)
from headflow import *

from aneurysm import Aneurysm

class PostProcessor(PostProcessorBase):
    def __init__(self):
        PostProcessorBase.__init__(self)
    



problem = Aneurysm()

scheme_pd = ParamDict(
                u_degree=1,
            )

scheme = IPCS(scheme_pd)
#scheme = SegregatedIPCS(scheme_pd)



postprocessor = PostProcessor()

wss = WSS(params=ParamDict(
            saveparams=ParamDict(
                save=True,
            )
        )
     )

postprocessor.add_field(wss)

pd = ParamDict(
    plot_solution=True,
)

solver = NSSolver(problem, scheme, postprocessor, params=pd)
solver.solve()