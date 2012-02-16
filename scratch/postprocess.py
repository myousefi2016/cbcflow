from __future__ import division
import numpy
import scipy
import pylab
from scipy import *
import os, re, sys, glob, math
from dolfin import *

number_of_steps = 3
casename = sys.argv[1]
results_dir = "results/success"
data_dir = "data/challenge"
case_dir = os.path.join(results_dir, casename)
assert os.path.exists(results_dir)
assert os.path.exists(data_dir)
assert os.path.exists(case_dir)
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
    # NB! Currently yielding the same function object each time!

    V = VectorFunctionSpace(mesh, "CG", 1)
    Vi = FunctionSpace(mesh, "CG", 1)
    u = Function(V)
    ui = Function(Vi)

    globpatterns = [os.path.join(case_dir, "u%d_at_t*.xml.gz" % i) for i in range(3)]
    regexps = [re.compile("u%d_at_t([0-9]*)_(.*).xml.gz" % i) for i in range(3)]

    filenames = {}
    times = {}
    for i in range(3):
        for fn in glob.glob(globpatterns[i]):
            m = regexps[i].search(fn)
            assert m
            ts, t = m.groups()
            ts = int(ts)
            filenames[ts] = filenames.get(ts,[]) + [fn]
            times[ts] = float(t)

    # Pick last few steps
    filenames = sorted(filenames.items())
    filenames = filenames[-number_of_steps:]
    print '\n'.join(map(str,filenames))

    n = ui.vector().size()
    for ts, fns in filenames:
        for i in range(3):
            f = File(fns[i])
            f >> ui.vector()
            u.vector()[i*n:(i+1)*n] = ui.vector() # NB! Assuming simple component layout!
        yield u, ts, times[ts]

def iterate_pressure_functions(mesh):
    # NB! Currently yielding the same function object each time!

    V = FunctionSpace(mesh, "CG", 1)
    p = Function(V)

    globpattern = os.path.join(case_dir, "p_at_t*.xml.gz")
    regexp = re.compile("p_at_t([0-9]*)_(.*).xml.gz")

    filenames = []
    times = {}
    for fn in glob.glob(globpattern):
        m = regexp.search(fn)
        assert m
        ts, t = m.groups()
        ts = int(ts)
        filenames.append((ts, fn))
        times[ts] = float(t)

    # Pick last few steps
    filenames = sorted(filenames)
    filenames = filenames[-number_of_steps:]
    print '\n'.join(map(str,filenames))

    for ts, fn in filenames:
        f = File(fn)
        f >> p.vector()
        yield p, ts, times[ts]

def evaluate_pressure_probes(cl, p):
    # Sample pressure at probe points from challenge readme1a
    probevalues = numpy.zeros((cl.shape[0],))
    for i in range(cl.shape[0]):
        probevalues[i] = p(numpy.array(cl[i,:3]))
    return probevalues

def postprocess(dowrite=False, doplot=False, dowritepvd=False):
    mesh = load_mesh()
    print mesh.num_vertices(), mesh.num_cells()

    cl = load_probes()
    cln = cl.shape[0]
    assert cl.shape[1] == 4
    x = numpy.array(cl[:,3])
    print "cl", cl.shape

    if doplot:
        line = None
        pylab.ion()

    if dowrite:
        if os.path.exists(probe_dir):
            print "NB! Already exists:", probe_dir
        else:
            print "Creating:", probe_dir
            os.mkdir(probe_dir)

    if dowritepvd:
        pvdfilename = os.path.join(probe_dir, "p.pvd")
        pvdfile = File(pvdfilename)
        upvdfilename = os.path.join(probe_dir, "u.pvd")
        upvdfile = File(upvdfilename)

    for p, ts, t in iterate_pressure_functions(mesh):
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
            probefilename = os.path.join(probe_dir, "p_t%.5e" % t)
            with open(probefilename, "w") as f:
                f.write('\n'.join(map(str,y)))
        if dowritepvd:
            pvdfile << p

    for u, ts, t in iterate_velocity_functions(mesh):
        if dowritepvd:
            upvdfile << u

    if doplot:
        pylab.ioff()
        pylab.show()

if __name__ == '__main__':
    args = sys.argv[2:]
    postprocess(doplot='p' in args,
                dowrite='w' in args,
                dowritepvd='v' in args)

