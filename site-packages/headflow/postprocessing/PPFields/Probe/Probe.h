#ifndef __PROBE_H
#define __PROBE_H

#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>

namespace dolfin
{
  class Probe
  {
  public:
      
    ~Probe();
    Probe(const Array<double>& x, const FunctionSpace& V);
    void eval(const Function& u);
    std::vector<double> get_probe(std::size_t i);
    std::size_t value_size();
    std::size_t number_of_evaluations();
    std::vector<double> coordinates();
    void erase(std::size_t i);
    void clear();
    
  private:
      
    std::vector<std::vector<double> > basis_matrix;
    std::vector<double> coefficients;
    double _x[3];
    boost::shared_ptr<const FiniteElement> _element;
    Cell* dolfin_cell;
    UFCCell* ufc_cell;
    std::size_t value_size_loc;    
    std::vector<std::vector<double> > _probes;
  };
}

#endif
