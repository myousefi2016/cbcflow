from __future__ import division
import numpy
import scipy
import pylab
from scipy import *
import os, re, sys, glob, math
from dolfin import *

casename = sys.argv[1]
results_dir = "results/success"
data_dir = "data/challenge"
case_dir = os.path.join(results_dir, casename)
probe_dir = os.path.join(case_dir, "probes")

def load_mesh(refinement=0, stationary=True, boundary_layers=True):
    if boundary_layers:
        mesh_filename = {
            0: "mesh_750k_BL_t.xml.gz",
            1: "mesh_2mio_BL_t.xml.gz",
            2: "mesh_4mio_BL_t.xml.gz",
            }[refinement]
    else:
        mesh_filename = {
            0: "mesh_500k.xml.gz",
            1: "mesh_1mio.xml.gz",
            2: "mesh_2mio.xml.gz",
            3: "mesh_4mio.xml.gz",
           }[refinement]
    mesh = Mesh(os.path.join(data_dir, mesh_filename))
    return mesh

def load_probes():
    cl = numpy.loadtxt(os.path.join(data_dir, "cl.dat"))
    return cl

def iterate_velocity_functions(mesh):
    # TODO: Invert this like for pressure, by using glob and re to find timesteps
    if 0:
        u = [u0,u1,u2]
        for i, ui in enumerate(u):
            fn = os.path.join(case_dir, "u%d_at_t_%.5e.xml" % (i,t))
            file = File(fn)
            file >> ui.vector()
    t = 0
    for u in []:
        yield u, t

def iterate_pressure_functions(mesh):
    # NB! Currently yielding the same function object each time!

    V = FunctionSpace(mesh, "CG", 1)
    p = Function(V)

    globpattern = os.path.join(case_dir, "p_at_t_*.xml")
    regexp = re.compile("p_at_t_(.*).xml")

    filenames = []
    for fn in glob.glob(globpattern):
        m = regexp.search(fn)
        assert m
        ts, = m.groups()
        t = float(ts)
        filenames.append((t, fn))

    # Pick last few steps
    filenames = sorted(filenames)
    filenames = filenames[-20:]
    print '\n'.join(map(str,filenames))

    for t, fn in filenames:
        f = File(fn)
        f >> p.vector()
        yield p, t

def evaluate_pressure_probes(cl, p):
    # Sample pressure at probe points from challenge readme1a
    probevalues = numpy.zeros((cl.shape[0],))
    for i in range(cl.shape[0]):
        probevalues[i] = p(numpy.array(cl[i,:3]))
    return probevalues

def postprocess(dowrite=False, doplot=False):
    mesh = load_mesh()
    print mesh.num_vertices(), mesh.num_cells()

    cl = load_probes()
    cln = cl.shape[0]
    assert cl.shape[1] == 4
    x = numpy.array(cl[:,3])
    print "cl", cl.shape

    line = None

    if dowrite:
        if os.path.exists(probe_dir):
            print "NB! Already exists:", probe_dir
        else:
            print "Creating:", probe_dir
            os.mkdir(probe_dir)

    for p, t in iterate_pressure_functions(mesh):
        y = evaluate_pressure_probes(cl, p)
        print "At t = %g, max dp = %g, probe dp = %g" % (t, p.vector().max()-p.vector().min(), max(y)-min(y))

        if doplot:
            if line is None:
                pylab.ion()
                line, = pylab.plot(x, y)
            else:
                line.set_data(y)
                pylab.draw()

        if dowrite:
            probefilename = os.path.join(probe_dir, "p_t%g" % t)
            with open(probefilename, "w") as f:
                f.write('\n'.join(map(str,self.probevalues)))

if __name__ == '__main__':
    postprocess(doplot=False, dowrite=False)
