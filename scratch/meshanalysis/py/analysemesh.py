#!/usr/bin/env python

import io
import numpy as np
import dolfin as df

from checklib import check, check_eq, check_lt, check_gt, check_le, check_ge, check_failures

from meshextractlib import extract_vertices, extract_domains, extract_domain_ids, replace_undefined_marker, count_marked_entities, count_subdomain_islands

def analyse_mesh_file(filename):
    # --- Load mesh file
    if filename == "UnitSquare":
        mesh = df.UnitSquare(2,1)
    else:
        mesh = df.Mesh(filename)

    # --- Order and initialize mesh entities
    ordered = mesh.ordered()
    mesh.init()
    ordered_after = mesh.ordered()
    check(ordered_after)

    # --- Get cell type # Still homogeneous for all dolfin meshes
    cell = mesh.ufl_cell()
    celltype = mesh.type()
    check_eq(type(celltype), df.CellType)
    check_eq(df.CellType.type2string(celltype.cell_type()), cell.domain())

    # --- Get geometry
    geometry = mesh.geometry()
    check_eq(type(geometry), df.MeshGeometry)

    coordinates = mesh.coordinates()
    check_eq(type(coordinates), np.ndarray)

    # --- Get topology
    topology = mesh.topology()
    check_eq(type(topology), df.MeshTopology)

    cells = mesh.cells()
    check_eq(type(cells), np.ndarray)

    # --- Get mesh data
    data = mesh.data()
    check_eq(type(data), df.MeshData)
    cell_domains, facet_domains = extract_domains(mesh)

    facet_domain_ids = extract_domain_ids(facet_domains)
    cell_domain_ids = extract_domain_ids(cell_domains)

    facet_undefined_id = replace_undefined_marker(facet_domains, facet_domain_ids)
    cell_undefined_id = replace_undefined_marker(cell_domains, cell_domain_ids)

    cell_domain_islands, num_cell_domain_islands = count_subdomain_islands(topology, cell_domains)
    facet_domain_islands, num_facet_domain_islands = count_subdomain_islands(topology, facet_domains)

    # --- Plot mesh and subdomains if wanted
    do_plot = 0
    if do_plot:
        df.plot(mesh, title='mesh')
        if 0 and facet_domains: # FIXME: Missing functionality in dolfin
            df.plot(facet_domains, title='facet_domains')
        if cell_domains:
            df.plot(cell_domains, title='cell_domains')
        df.interactive()

    # --- Get dimensions
    top_dim = topology.dim()
    geom_dim = geometry.dim()
    check_eq(top_dim, geom_dim) # Will not be true for manifolds

    num_cells = mesh.num_cells()
    num_facets = mesh.num_facets()
    num_vertices = mesh.num_vertices()

    check_eq(num_vertices, mesh.size(0))
    check_eq(num_facets, mesh.size(top_dim-1))
    check_eq(num_cells, mesh.size(top_dim))
    check_eq(coordinates.shape, (num_vertices, geom_dim))

    # --- Compute subdomain related geometric quantities
    one = df.Constant(1.0)

    interior_areas = [df.assemble(df.avg(one)*df.dS(d), mesh=mesh) for d in facet_domain_ids]
    # Not implemented in FFC:
    #cell_areas = [df.assemble(cell.surface_area/cell.volume*df.dx(d), mesh=mesh)
    #              for d in cell_domain_ids]

    # --- Compute coordinate related geometric quantities
    bboxmin = tuple(min(coordinates[:,j]) for j in range(geom_dim))
    bboxmax = tuple(max(coordinates[:,j]) for j in range(geom_dim))
    bboxvol = df.product(bboxmax[j]-bboxmin[j] for j in range(geom_dim))
    # TODO: This barycenter is inaccurate for graded meshes:
    barycenter = tuple(sum(coordinates) / len(coordinates))
    bary_coords = coordinates - barycenter
    bary_radius = tuple(max(abs(bary_coords[:,j])) for j in range(geom_dim))
    bary_bboxmin = tuple(min(bary_coords[:,j]) for j in range(geom_dim))
    bary_bboxmax = tuple(max(bary_coords[:,j]) for j in range(geom_dim))

    # --- Build barycenter, bbox, radius, etc. for each cell subdomain
    cell_domain_data = {}
    for dind, d in enumerate(cell_domain_ids):
        # Integrate 1 and x over subdomain to find volume and barycenter
        sd_volume = df.assemble(one*df.dx(d), mesh=mesh)
        sd_volume = sd_volume or 1 # FIXME: Better way to handle undefined domains?
        sd_unscaled_barycenter = tuple(df.assemble(cell.x[j]*df.dx(d), mesh=mesh)
                                       for j in range(geom_dim))
        sd_barycenter = tuple(sd_unscaled_barycenter[j] / sd_volume
                              for j in range(geom_dim))

        # Find coordinates in this subdomain
        sd_vertices = extract_vertices(topology, cell_domains, d)
        sd_coords = coordinates[sd_vertices,:]

        # Build bbox from global coordinates
        sd_bboxmin = tuple(min(sd_coords[:,j]) for j in range(geom_dim))
        sd_bboxmax = tuple(max(sd_coords[:,j]) for j in range(geom_dim))

        # Build bboxmin, bboxmax, radius from barycentric relative coordinates
        sd_bary_coords = sd_coords - sd_barycenter
        sd_bary_radius = max(np.sqrt(sum(sd_bary_coords[i,:]**2))
                             for i in xrange(sd_bary_coords.shape[0]))
        sd_bary_bboxmin = tuple(min(sd_bary_coords[:,j]) for j in range(geom_dim))
        sd_bary_bboxmax = tuple(max(sd_bary_coords[:,j]) for j in range(geom_dim))

        # Store data for later printing
        cell_domain_data[d] = { 'bary_radius': sd_bary_radius,
                                'bary_bboxmin': sd_bary_bboxmin,
                                'bary_bboxmax': sd_bary_bboxmax,
                                'bboxmin': sd_bboxmin,
                                'bboxmax': sd_bboxmax,
                                'barycenter': sd_barycenter,
                                }

    # --- Build barycenter, bbox, radius, etc. for each facet subdomain
    facet_domain_data = {}
    for dind, d in enumerate(facet_domain_ids):
        # Integrate 1 and x over subdomain to find area and barycenter
        sd_area = df.assemble(one*df.ds(d), mesh=mesh)
        sd_area = sd_area or 1
        sd_unscaled_barycenter = tuple(df.assemble(cell.x[j]*df.ds(d), mesh=mesh)
                                       for j in range(geom_dim))
        sd_barycenter = tuple(sd_unscaled_barycenter[j] / sd_area
                              for j in range(geom_dim))

        # Find coordinates in this subdomain
        sd_vertices = extract_vertices(topology, facet_domains, d)
        sd_coords = coordinates[sd_vertices,:]

        # Build bbox from global coordinates
        sd_bboxmin = tuple(min(sd_coords[:,j]) for j in range(geom_dim))
        sd_bboxmax = tuple(max(sd_coords[:,j]) for j in range(geom_dim))

        # Build bboxmin, bboxmax, radius from barycentric relative coordinates
        sd_bary_coords = sd_coords - sd_barycenter
        sd_bary_radius = max(np.sqrt(sum(sd_bary_coords[i,:]**2))
                             for i in xrange(sd_bary_coords.shape[0]))
        sd_bary_bboxmin = tuple(min(sd_bary_coords[:,j]) for j in range(geom_dim))
        sd_bary_bboxmax = tuple(max(sd_bary_coords[:,j]) for j in range(geom_dim))

        # Store data for later printing
        facet_domain_data[d] = { 'bary_radius': sd_bary_radius,
                                 'bary_bboxmin': sd_bary_bboxmin,
                                 'bary_bboxmax': sd_bary_bboxmax,
                                 'bboxmin': sd_bboxmin,
                                 'bboxmax': sd_bboxmax,
                                 'barycenter': sd_barycenter,
                                 }

    # --- Compute cell related geometric quantities
    # TODO: For each cell, compute:
    # - volume
    # - facet areas
    # - total surface area
    # - smallest surface area
    # - largest surface area
    # - average surface area
    # - circumradius
    # - smallest angle
    # - largest angle
    # TODO: Build statistics of these quantities

    hmin = mesh.hmin()
    hmax = mesh.hmax()
    check_le(hmin, hmax)
    check_gt(hmin, 0.0)

    # --- Consistency computations with subdomains

    cell_marker_counts = count_marked_entities(cell_domains)
    facet_marker_counts = count_marked_entities(facet_domains)

    counted_cells = 0
    for d in cell_domain_ids:
        # The integral of 1/|K| over all cells K should be the number of cells
        count = int(df.assemble(1/cell.volume*df.dx(d), mesh=mesh))
        check_eq(count, cell_marker_counts[d])
        counted_cells += count

    counted_exterior_facets = 0
    counted_interior_facets = 0
    for d in facet_domain_ids:
        # The integral of 1/|T| over all facets T should be the number of facets
        exterior_count = int(df.assemble(one/cell.facet_area*df.ds(d), mesh=mesh))
        interior_count = int(df.assemble(df.avg(one/cell.facet_area)*df.dS(d), mesh=mesh))
        check_eq(exterior_count + interior_count, facet_marker_counts[d])
        print 'dei', d, exterior_count, interior_count
        check(exterior_count == 0 or interior_count == 0)
        counted_exterior_facets += exterior_count
        counted_interior_facets += interior_count

    # TODO: Get number of cells in each subdomain and compare separately
    check_eq(counted_cells, num_cells)
    # TODO: Get number of facets in each subdomain and compare separately
    check_eq(counted_interior_facets+counted_exterior_facets, num_facets)
    facets_per_cell = top_dim + 1 # TODO: simplex specific!
    check_eq(2*counted_interior_facets+counted_exterior_facets, facets_per_cell*num_cells)

    # For each set of interior facets, the normal components in + and - directions
    # should sum up to the opposite number
    for d in facet_domain_ids:
        for c in range(geom_dim):
            interior_normals = [df.assemble(cell.n(r)[c]*df.dS(d), mesh=mesh)
                                for r in ('+', '-')]
            check_lt(abs(interior_normals[0] + interior_normals[1]), 1e-10)

    # For each geometric dimension, the normal component should integrate
    # to 0 over the entire boundary
    for c in range(geom_dim):
        boundary_normals = [df.assemble(cell.n[c]*df.ds(d), mesh=mesh)
                            for d in facet_domain_ids]
        print boundary_normals
        print sum(boundary_normals)
        check_lt(abs(sum(boundary_normals)), 1e-9)

    # --- Print summary
    mesh_vars = [
        ('Mesh',
         ('filename',
          'ordered',)),
        ('Dimensions',
         ('num_vertices', 'num_cells',  'num_facets',
          'counted_cells', 'counted_interior_facets', 'counted_exterior_facets'
          )),
        ('Coordinates',
         ('bboxmin', 'bboxmax', 'bboxvol',
          'barycenter', 'bary_radius', 'bary_bboxmin', 'bary_bboxmax',)),
        ('Cell geometry',
         ('hmin', 'hmax',)),
        ('Domain geometry',
         ('interior_areas',)),
        ('Domain markers',
         ('facet_undefined_id', 'cell_undefined_id',
#          'cell_domain_islands', 'facet_domain_islands',
          'num_cell_domain_islands', 'num_facet_domain_islands',)),
        ]
    print_variables(mesh_vars, locals())

    subdomain_varnames = list(cell_domain_data[iter(cell_domain_ids).next()].keys())
    for d in cell_domain_ids:
        subdomain_vars = [("Cell subdomain %d" % d, subdomain_varnames)]
        print_variables(subdomain_vars, cell_domain_data[d])

    subdomain_varnames = list(facet_domain_data[iter(facet_domain_ids).next()].keys())
    for d in facet_domain_ids:
        subdomain_vars = [("Facet subdomain %d" % d, subdomain_varnames)]
        print_variables(subdomain_vars, facet_domain_data[d])

    failures = check_failures()
    if failures:
        print_failures(failures)

def print_variables(printme, loc):
    for caption, keys in printme:
        print "\n--- %s" % (caption,)
        n = max(len(k) for k in keys)
        for k in keys:
            sp = " "*(n-len(k))
            print "%s%s = %s" % (k, sp, loc[k])

def print_failures(failures):
    print
    print "--- Some checks failed, summarized here:"
    for f in failures:
        print
        print f
    print

def main(args):
    for fn in args:
        analyse_mesh_file(fn)
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv[1:]))
