#ifndef __UNASSEMBLED_H
#define __UNASSEMBLED_H

#include <dolfin/fem/UFC.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/fem/Form.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/fem/Assembler.h>
#include <petscmat.h>
#include <string>

//#include <typeinfo>
//#include <stdexcept>

namespace dolfin
{
    class UnassembledMatrix
    {
    public:
        UnassembledMatrix(const Form& a);
        //UnassembledMatrix(const Form& a, std::shared_ptr<Function> coeff);
        //UnassembledMatrix(const Form& a, std::vector<std::shared_ptr<Function> > coeffs);
        virtual ~UnassembledMatrix();
        //virtual void clear();
        //dolfin::Matrix& assemble();
        //void assemble(dolfin::GenericTensor& A, dolfin::Function& coeff);
        void assemble(dolfin::GenericTensor& A);
        void assemble(dolfin::GenericTensor& A, std::vector<std::shared_ptr<const Function> > coeffs);
    
    protected:
        //void init(std::vector<std::size_t> coeff_indices);
        void init();
        std::vector<std::size_t> get_function_indices();
        
        std::vector<std::vector<dolfin::la_index> > indices;
        std::vector<double> all_values;
        

        const Form _a;
        const std::vector<std::size_t> coeff_indices;
        const std::size_t _num_coefficients;
        
        std::shared_ptr<const GenericDofMap> dm;
        std::size_t _max_cols;
        
    };
}

#endif