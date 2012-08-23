
from itertools import chain
import numpy as np
from checklib import check, check_eq, check_lt

def extract_vertices(topology, markers, value):
    """
    Return array of vertex indices in topology covered by markers with given value.
    """
    # Probably faster for large subdomains to use arrays
    # instead of set, to get O(n) instead of O(n log n)
    conn = topology(markers.dim(), 0)
    ma = markers.array()
    mv, = np.where(ma == value)
    vertices = set(chain(*[conn(i) for i in mv]))
    return np.asarray(list(vertices))

def extract_domains(mesh):
    """
    domains: MeshDomains
    """
    td = mesh.topology().dim()
    domains = mesh.domains()
    check_eq(type(domains), df.MeshDomains)

    cell_domains = domains.cell_domains(mesh)
    if cell_domains is None:
        cell_domains = df.MeshFunction('uint', mesh, td)
        cell_domains.set_all(0)
    check_eq(type(cell_domains), df.MeshFunctionUInt)

    facet_domains = domains.facet_domains(mesh)
    if facet_domains is None:
        facet_domains = df.MeshFunction('uint', mesh, td-1)
        facet_domains.set_all(0)
    check_eq(type(facet_domains), df.MeshFunctionUInt)

    return cell_domains, facet_domains

def extract_domain_ids(meshfunction):
    # Cumbersome way to get set of markers, easier way?
    if meshfunction:
        domain_ids = set(meshfunction.array())
    else:
        domain_ids = set((0,))
    return domain_ids

def replace_undefined_marker(markers, ids):
    undef = df.MeshDomains.default_unset_value
    if undef in ids and len(ids) > 1:
        ids.remove(undef)
        largest = max(ids)
        newundef = largest+1
        ids.add(newundef)
        ma = markers.array()
        ma[np.where(ma == undef)] = newundef
        return newundef
    else:
        return undef

def count_marked_entities(markers):
    counts = {}
    for i in markers.array():
        counts[i] = counts.get(i,0) + 1
    return counts

def count_subdomain_islands(topology, markers):
    td = topology.dim()
    md = markers.dim()
    n = markers.size()

    conn = topology(md, md)

    islands = []
    visited = np.zeros((n,), dtype=np.intc)

    all_visited_below = 0
    num_visited = 0
    while num_visited < n:
        for i in xrange(all_visited_below, n):
            if not visited[i]:
                break
        all_visited_below = i
        todo = [all_visited_below]
        current_island = []

        while todo:
            i = todo.pop()
            if visited[i]:
                continue

            visited[i] = 1
            num_visited += 1
            current_island.append(i)

            mi = markers[i]
            c = conn(i)
            for j in c:
                if not visited[j] and mi == markers[j]:
                    todo.append(j)

        islands.append((mi, current_island))

    # Build marker id -> list of islands
    subdomain_islands = {}
    for mi, isl in islands:
        subdomain_islands[mi] = subdomain_islands.get(mi,[]) + [isl]

    # Build marker id -> number of islands
    num_subdomain_islands = {}
    for mi, isls in subdomain_islands.iteritems():
        num_subdomain_islands[mi] = len(isls)

    return subdomain_islands, num_subdomain_islands

import dolfin as df
import unittest
class TestMeshAnalysis(unittest.TestCase):
    def test_extract_vertices(self):
        mesh = df.UnitSquare(1,1)
        mesh.init()
        top = mesh.topology()
        geom = mesh.geometry()
        coordinates = mesh.coordinates()

        resetvalue = 0
        markervalue = 1

        for c in range(mesh.num_cells()):
            print "cell", c
            cmf = df.MeshFunction('uint', mesh, top.dim())
            cmf.set_all(resetvalue)
            cmf[c] = markervalue
            cv = extract_vertices(top, cmf, markervalue)
            print cv
            coords = coordinates[cv,:]
            print coords

        for f in range(mesh.num_facets()):
            print "facet", f
            fmf = df.MeshFunction('uint', mesh, top.dim()-1)
            fmf.set_all(resetvalue)
            fmf[f] = markervalue
            fv = extract_vertices(top, fmf, markervalue)
            print fv
            coords = coordinates[fv,:]
            print coords

        if 0:
            df.plot(mesh)
            df.interactive()

if __name__ == "__main__":
    unittest.main()
