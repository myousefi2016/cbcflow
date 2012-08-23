// Copyright (C) 2006 Ola Skavhaug
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// Modified by Anders Logg, 2007.
//
// First added:  2006-11-29
// Last changed: 2012-07-05

#include <dolfin.h>
#include <cmath>
#include <cassert>
#include <limits>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>

#include "checklib.h"
#include "loglib.h"

template<typename T>
class Accumulator
{
public:
    Accumulator():
        min(std::numeric_limits<T>::max()),
        max(std::numeric_limits<T>::min()),
        avg(0),
        sum(0),
        num(0)
    {
    }

    void update(const T & value)
    {
        ++num;
        sum += value;
        min = std::min(min, value);
        max = std::max(max, value);
    }

    void finalize()
    {
        avg = sum / T(num);
    }

    T min;
    T max;
    T avg;
    T sum;
    int num;
};

template<typename T>
class Histogram
{
public:
    Histogram(int bins, T min, T max):
        counts(bins), percents(bins), min(min), max(max), total(0), outside(0)
    {}

    void update(const T & value)
    {
        int hist_ind = int(counts.size() * (value-min) / (max-min));
        if (hist_ind < 0 || hist_ind >= counts.size())
        {
            std::cerr << "ERROR "
                      << "hist_ind " << hist_ind
                      << ", value " << value
                      << ", min " << min
                      << ", max " << max << std::endl;
            ++outside;
        }
        else
        {
            counts[hist_ind] += 1;
        }
        ++total;
    }

    void finalize()
    {
        for (int i = 0; i < percents.size(); ++i)
            percents[i] = float(counts[i]) / float(total);
    }

    std::vector<int> counts;
    std::vector<float> percents;
    const T min;
    const T max;
    int total;
    int outside;
};


template<typename T>
class Histogram2
{
public:
    Histogram2(int bins):
        counts(bins),
        percents(bins),
        min(std::numeric_limits<T>::max()),
        max(std::numeric_limits<T>::min()),
        total(0)
    {}

    void update(const T & value)
    {
        values.push_back(value);
        min = std::min(min, value);
        max = std::max(max, value);
        ++total;
    }

    std::pair<T,T> range(int bin)
    {
        T from = min + (max-min) * bin / counts.size();
        T to = min + (max-min) * (bin+1) / counts.size();
        return std::make_pair(from, to);
    }

    void finalize()
    {
        for (typename values_container::iterator iter = values.begin(); iter != values.end(); ++iter)
        {
            T value = *iter;
            assert(value >= min);
            assert(value <= max);
            int hist_ind = int(counts.size() * (value-min) / (max-min));
            if (value == max && hist_ind == counts.size())
                hist_ind--;
            assert(hist_ind >= 0);
            assert(hist_ind < counts.size());
            counts[hist_ind] += 1;
        }
        for (int i = 0; i < percents.size(); ++i)
            percents[i] = float(counts[i]) / float(total);
    }

    typedef std::list<T> values_container;
    values_container values;
    std::vector<int> counts;
    std::vector<float> percents;
    T min;
    T max;
    int total;
};



void analyse_mesh(dolfin::Mesh & mesh)
{
    dolfin::MeshGeometry & geometry = mesh.geometry();
    dolfin::MeshTopology & topology = mesh.topology();
    dolfin::MeshData & data = mesh.data();
    dolfin::CellType & celltype = mesh.type();
    bool initially_ordered = mesh.ordered();
    std::string cell_type_name = dolfin::CellType::type2string(celltype.cell_type());
    std::string facet_type_name = dolfin::CellType::type2string(celltype.facet_type());

    // Check ordering
    log2("initially_ordered", initially_ordered);
    mesh.init();
    check(mesh.ordered());

    int td = topology.dim();
    int gd = geometry.dim();

    log2("Stats", "Dimensions");
    log2("top dim", td);
    log2("geom dim", gd);
    log2("cell", cell_type_name);
    log2("facet", facet_type_name);
    log2("num vertices", topology.size(0));
    log2("num facets", topology.size(td-1));
    log2("num cells", topology.size(td));

    // Check for dimension consistency (not likely to fail, written mostly to get overview of interface!)
    check_eq(topology.dim(), celltype.dim());
    /*
    switch (topology.dim())
    {
    case 1:
        check_eq(celltype.cell_type(), dolfin::CellType::interval);
        check_eq(cell_type_name, "interval");
        check_eq(facet_type_name, "vertex");
        break;
    case 2:
        check_eq(celltype.cell_type(), dolfin::CellType::triangle);
        check_eq(cell_type_name, "triangle");
        check_eq(facet_type_name, "interval");
        break;
    case 3:
        check_eq(celltype.cell_type(), dolfin::CellType::tetrahedron);
        check_eq(cell_type_name, "tetrahedron");
        check_eq(facet_type_name, "triangle");
        break;
    }
    */

    const double pi = acos(-1.0);

    // Initialize statistics
    Accumulator<double> area_stats;
    Histogram2<double> area_hist(10);
    Accumulator<double> edge_stats;
    Histogram2<double> edge_hist(10);
    Accumulator<double> angle_stats;
    Histogram<double> angle_hist(10, 0.0, pi);

    double sum_nn = 0.0;
    int num_nn = 0;

    dolfin::CellFunction<double> cf_volume(mesh);
    dolfin::CellFunction<double> cf_orientation(mesh);
    dolfin::CellFunction<double> cf_diameter(mesh);
    dolfin::CellFunction<double> cf_max_area_ratio(mesh);
    dolfin::CellFunction<double> cf_max_edge_ratio(mesh);
    dolfin::CellFunction<double> cf_max_midangle(mesh);

    // Iterate over all cells
    for (dolfin::CellIterator cell(mesh); !cell.end(); ++cell)
    {
        dolfin::Cell & c = *cell;

        // Get cell indexing quantities
        int ind = c.index();
        int dim = c.dim();
        int num_facets = c.num_entities(dim-1);
        int num_edges = c.num_entities(1);
        int num_vertices = c.num_entities(0);

        // Compute cell quantities
        double volume = c.volume();
        double orientation = c.orientation();
        double diameter = c.diameter();

        // Extrema of per-facet or per-edge quantities
        double max_area_ratio = 0.0;
        double max_edge_ratio = 0.0;
        double max_midangle = 0.0;

        // Compute per-facet quantities
        for (int facet = 0; facet < num_facets; ++facet)
        {
            double area = c.facet_area(facet);
            dolfin::Point n = c.normal(facet);
            double nn = n.norm();

            sum_nn += nn;
            num_nn++;

            // Compute facet-facet quantities
            for (int facet2 = 1; facet2 < num_facets; ++facet2)
            {
                if (facet < facet2)
                {
                    double area2 = c.facet_area(facet2);
                    double area_ratio = area2 / area;
                    if (area_ratio < 1.0)
                        area_ratio = 1.0 / area_ratio;
                    area_stats.update(area_ratio);
                    area_hist.update(area_ratio);
                    max_area_ratio = std::max(max_area_ratio, area_ratio);

                    dolfin::Point n2 = c.normal(facet2);
                    double nn2 = n.norm();
                    double arcangle = n.dot(n2) / (nn2*nn);
                    double angle = acos(arcangle);
                    double midangle = std::abs(angle-pi/2.0) * 360.0/pi;
                    angle_stats.update(angle);
                    angle_hist.update(angle);
                    max_midangle = std::max(max_midangle, midangle);
                }
            }
        }

        if (dim > 2)
        {
            // Compute per-edge quantities
            const dolfin::uint * vertices = c.entities(0);
            for (int vertex = 0; vertex < num_vertices; ++vertex)
            {
                for (int vertex2 = 1; vertex2 < num_vertices; ++vertex2)
                {
                    if (vertex < vertex2)
                    {
                        dolfin::Point p1(dim, geometry.x(vertices[vertex]));
                        dolfin::Point p2(dim, geometry.x(vertices[vertex2]));
                        dolfin::Point e = p2 - p1;
                        double en = e.norm();
                        double edge_ratio = en / diameter;
                        edge_stats.update(edge_ratio);
                        edge_hist.update(edge_ratio);
                        max_edge_ratio = std::max(max_edge_ratio, edge_ratio);
                    }
                }
            }
        }

        // Store results per cell
        cf_volume[ind] = volume;
        cf_orientation[ind] = orientation;
        cf_diameter[ind] = diameter;
        cf_max_area_ratio[ind] = max_area_ratio;
        cf_max_edge_ratio[ind] = max_edge_ratio;
        cf_max_midangle[ind] = max_midangle;
    }
    // Finalize statistics
    area_stats.finalize();
    area_hist.finalize();
    angle_stats.finalize();
    angle_hist.finalize();
    edge_stats.finalize();
    edge_hist.finalize();

    double avg_nn = sum_nn / num_nn;

    // Store to file
    dolfin::File cf_volume_file("cf_volume.pvd");
    cf_volume_file << cf_volume;
    dolfin::File cf_orientation_file("cf_orientation.pvd");
    cf_orientation_file << cf_orientation;
    dolfin::File cf_diameter_file("cf_diameter.pvd");
    cf_diameter_file << cf_diameter;
    dolfin::File cf_max_area_ratio_file("cf_max_area_ratio.pvd");
    cf_max_area_ratio_file << cf_max_area_ratio;
    dolfin::File cf_max_edge_ratio_file("cf_max_edge_ratio.pvd");
    cf_max_edge_ratio_file << cf_max_edge_ratio;
    dolfin::File cf_max_midangle_file("cf_max_midangle.pvd");
    cf_max_midangle_file << cf_max_midangle;

    // Print to screen
    log2("Stats: ", "nn");
    log2("  num: ", num_nn);
    log2("  avg: ", avg_nn);

    log2("Stats: ", "Area ratio");
    log2("  num: ", area_stats.num);
    log2("  min: ", area_stats.min);
    log2("  max: ", area_stats.max);
    log2("  avg: ", area_stats.avg);

    log2("Histogram", "Area ratio");
    for (int i=0; i<area_hist.counts.size(); ++i)
    {
        std::pair<double,double> r = area_hist.range(i);
        log5(i, r.first, r.second, area_hist.counts[i], area_hist.percents[i]);
    }

    log2("Stats: ", "Edge ratio");
    log2("  num: ", edge_stats.num);
    log2("  min: ", edge_stats.min);
    log2("  max: ", edge_stats.max);
    log2("  avg: ", edge_stats.avg);

    log2("Histogram", "Edge ratio");
    for (int i=0; i<edge_hist.counts.size(); ++i)
    {
        std::pair<double,double> r = edge_hist.range(i);
        log5(i, r.first, r.second, edge_hist.counts[i], edge_hist.percents[i]);
    }

    log2("Stats: ", "Angle");
    log2("  num: ", angle_stats.num);
    log2("  min: ", angle_stats.min);
    log2("  max: ", angle_stats.max);
    log2("  avg: ", angle_stats.avg);

    log2("Histogram", "Angle");
    for (int i=0; i<angle_hist.counts.size(); ++i)
        log3(i, angle_hist.counts[i], angle_hist.percents[i]);
}


void usage()
{
    std::cout << "Usage: analysemesh <meshfilename>" << std::endl;
}

int main(int argc, const char * argv[])
{
    std::vector<std::string> args(argv+1, argv+argc);

    if (args.size() != 1)
    {
        usage();
        return -1;
    }

    std::string mesh_filename(args[0]);

    dolfin::Mesh mesh(mesh_filename);

    analyse_mesh(mesh);

    return 0;
}
