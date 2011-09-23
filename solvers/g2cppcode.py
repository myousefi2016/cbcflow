__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-01-01"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

cppcode_d1 = """
class Delta1 : public Expression
{
public:

  Delta1() : Expression() {}

  //void eval(Array<double>& values, const Data& data) const
  void eval(Array<double>& values, const Array<double>& x,  const ufc::cell& cell) const
  {
    const uint D = cell.topological_dimension;
    const uint cell_index = cell.entity_indices[D][0];
    values[0] = _values[cell_index];
    //values[0] = _values[data.cell().index()];
  }

  void update(boost::shared_ptr<dolfin::Function> u, double nu, double dt, double C1)
  {
    const Mesh& mesh = *u->function_space().mesh();
    if (_values.size() != mesh.num_cells())
      _values.resize(mesh.num_cells());

    const GenericVector& x = u->vector();
    const uint gdim = mesh.geometry().dim();

    for (CellIterator cell(mesh); !cell.end(); ++cell)
    {
      double U = 0.0;
      for (VertexIterator vertex(*cell); !vertex.end(); ++vertex)
      {
        for (uint i = 0; i < gdim; ++i)
        {
          const double Ui = x[i*mesh.num_vertices() + vertex->index()];
          U += Ui*Ui;
        }
      }
      // FIXME: Should reverse order of sqrt and /=, but this is how it's done in Unicorn
      U = sqrt(U);
      U /= static_cast<double>(gdim + 1);

      const double h = cell->diameter();
      const uint i = cell->index();

      if (h/nu > 1.0 || nu < 1.0e-10)
        _values[i] = C1 * (0.5/sqrt(1.0/(dt*dt) + (U/h)*(U/h)));
      else
        _values[i] = C1*h*h;
    }
  }

private:

  std::vector<double> _values;

};"""

cppcode_d2 = """
class Delta2 : public Expression
{
public:

  Delta2() : Expression() {}

  //void eval(Array<double>& values, const Data& data) const
  void eval(Array<double>& values, const Array<double>& x,  const ufc::cell& cell) const
  {
    const uint D = cell.topological_dimension;
    const uint cell_index = cell.entity_indices[D][0];
    values[0] = _values[cell_index];
    //values[0] = _values[data.cell().index()];
  }

  void update(boost::shared_ptr<dolfin::Function> u, double nu, double dt, double C2)
  {
    const Mesh& mesh = *u->function_space().mesh();
    if (_values.size() != mesh.num_cells())
      _values.resize(mesh.num_cells());

    for (CellIterator cell(mesh); !cell.end(); ++cell)
    {
      const double h = cell->diameter();
      const uint i = cell->index();

      if (h/nu > 1.0 || nu < 1.0e-10)
        _values[i] = C2*h;
      else
        _values[i] = C2*h*h;
    }
  }

private:

  std::vector<double> _values;

};"""
