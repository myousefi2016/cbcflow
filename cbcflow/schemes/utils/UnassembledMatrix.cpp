#include "UnassembledMatrix.h"

using namespace dolfin;

UnassembledMatrix::~UnassembledMatrix()
{
}

//UnassembledMatrix::UnassembledMatrix(const Form& a) : _a(a), _num_coefficients(a.num_coefficients())
UnassembledMatrix::UnassembledMatrix(const Form& a) : _a(a), coeff_indices(get_function_indices()), _num_coefficients(coeff_indices.size())
{
    //std::vector<std::size_t> coeff_indices(_num_coefficients);
    //std::iota(coeff_indices.begin(), coeff_indices.end(), 0);
    //init(coeff_indices);
    std::cout << "Num coefficients (functions): " << _num_coefficients << std::endl;
    
    init();
}

std::vector<std::size_t> UnassembledMatrix::get_function_indices()
{
    std::vector<std::size_t> function_indices;
    for (std::size_t i=0; i<_a.num_coefficients(); i++)
    {
        std::shared_ptr<const dolfin::Function> func = std::dynamic_pointer_cast<const dolfin::Function>(_a.coefficients()[i]);
        if (func != NULL)
        {
            function_indices.push_back(i);
            std::cout << "Index: " << i << std::endl;
            std::cout << (func==NULL) << std::endl;
        }
        //coeffs[i] = std::dynamic_pointer_cast<const dolfin::Function>(_a.coefficients()[i]);
        //std::cout << (coeffs[i] == NULL) << std::endl;
    }
    //assemble(A, coeffs);
    return function_indices;
}

/*
UnassembledMatrix::UnassembledMatrix(const Form& a, std::shared_ptr<Function> coeff) : _a(a), _num_coefficients(1)
{
    std::cout << "In UnassembledMatrix constructor" << std::endl;
    
}
*/
/*
UnassembledMatrix::UnassembledMatrix(const Form& a, std::vector<std::shared_ptr<Function> > coeffs) : _a(a), _num_coefficients(coeffs.size())
{
    std::cout << "In UnassembledMatrix constructor" << std::endl;
    
}
*/


void UnassembledMatrix::init()
//void UnassembledMatrix::init(std::vector<std::size_t> coeff_indices)
//UnassembledMatrix::UnassembledMatrix(const Form& a, std::vector<std::size_t> coeff_indices) : _a(a)
{
    //std::cout << (std::string(DOLFIN_VERSION)==std::string("1.4.0")) << std::endl;
    //std::cout << std::string(DOLFIN_VERSION).compare(std::string("1.4.0")) << std::endl;
    //std::cout << string::compare(DOLFIN_VERSION, "1.4.0") << std::endl;
    
    _max_cols = 0;
    std::cout << "In UnassembledMatrix constructor" << std::endl;
    const Mesh& mesh = _a.mesh();
    const std::size_t D = mesh.topology().dim();
    std::cout << "Dimension: " << D << std::endl;
    UFC ufc(_a);
    
    std::shared_ptr<const FunctionSpace> V = _a.function_space(0);
    dm = V->dofmap();
    //std::size_t N_unowned = dm->local_to_global_unowned().size();

    // Get number of dofs on this process (owned+unowned)
    std::size_t N = dm->ownership_range().second - dm->ownership_range().first + dm->local_to_global_unowned().size(); // 1.4.0+
    // Create space for number of dofs
    indices.resize(N);
    /*
    std::size_t K = 0;
    std::cout << "Num owned dofs: " << dm->ownership_range().second - dm->ownership_range().first << std::endl;
    std::cout << "Num unowned dofs: " << dm->local_to_global_unowned().size() << std::endl;
    std::cout << "Total num local dofs: " << N << std::endl;
    std::cout << MPI::rank(MPI_COMM_WORLD) << std::endl;
    */

    std::vector< std::vector<Cell> > dof_cell_support(N);

    // Find map between dofs and the cells the dof basis function has support on
    for (CellIterator cell(mesh); !cell.end(); ++cell)
    {
        const std::vector<dolfin::la_index>& dofs = dm->cell_dofs(cell->index());

        for (std::size_t i = 0; i < dofs.size(); i++)
            dof_cell_support[dofs[i]].push_back(*cell);

    }

    // Create dof matrices by integrating over cells in dofs support
    ufc::cell ufc_cell;
    ufc::cell_integral* integral = ufc.default_cell_integral.get();
    
    // Form rank
    const std::size_t form_rank = ufc.form.rank();

      // Collect pointers to dof maps
    std::vector<const GenericDofMap*> dofmaps;
    for (std::size_t i = 0; i < form_rank; ++i)
        dofmaps.push_back(_a.function_space(i)->dofmap().get());

    // Vector to hold dof map for a cell
    std::vector<std::vector<dolfin::la_index>> dofs(form_rank);
    std::vector<dolfin::la_index> local_dofs;
    
    std::vector<double> zero_vector(coeff_indices.size());
    std::fill(zero_vector.begin(), zero_vector.end(), 0.0);
    
    //std::pair<std::size_t, std::size_t> ij;
    
    std::size_t row, column, coeff_index;
    std::vector<dolfin::la_index> columns;
    
    std::vector<double> vertex_coordinates(D);
    
    std::size_t tmp_idx;
    
    // Create container (and pointer) for dof values at each cell    
    std::vector<std::vector<double> > _w;
    _w.resize(_num_coefficients);
    std::vector<double*> w_pointer;
    w_pointer.resize(_num_coefficients);
    for (std::size_t i=0; i<_num_coefficients; i++)
    {
        _w[i].resize(dm->cell_dimension(0));
        w_pointer[i] = &_w[i][0];
    }

    // Local indices for local dense matrix
    std::size_t dof_row_index, dof_col_index, dof_flat_index;
    std::vector<double> dofmatrix_flat;
    
    // Container for holding dof indices in local cell support
    std::vector<dolfin::la_index> dof_indices;
    dof_indices.reserve(1000);
    
    // Map from function space dof to renumbered local dof
    std::map<dolfin::la_index, std::size_t> dofmatrix_to_local_map;
    
    // Index of dof in cell
    dolfin::la_index cell_dof_index;
    
    // Create vector to hold temporary values
    std::vector<double> tmp_values(coeff_indices.size());
    std::vector<std::vector<double> > values;
    
    
    for (dolfin::la_index dof=0; dof<N; dof++) {
        // Find all dofs in cell support of this dof
        dof_indices.clear();
        for(std::vector<Cell>::iterator cell = dof_cell_support[dof].begin(); cell != dof_cell_support[dof].end(); ++cell)
        {
            const std::vector<dolfin::la_index>& celldofs = dm->cell_dofs(cell->index());
            dof_indices.insert(dof_indices.end(), celldofs.begin(), celldofs.end());
        }
        std::sort(dof_indices.begin(), dof_indices.end()); 
        dof_indices.erase( std::unique( dof_indices.begin(), dof_indices.end() ), dof_indices.end() );

        // Create local numbering of dofs (0, 1, ..., dof_indices.size()-1)
        dofmatrix_to_local_map.clear();
        std::size_t local_idx=0;
        for (std::size_t j=0; j<dof_indices.size(); j++)
            dofmatrix_to_local_map[dof_indices[j]] = local_idx++;
        
        // Make space for dense local matrix
        dofmatrix_flat.clear();
        dofmatrix_flat.resize(dof_indices.size()*dof_indices.size()*coeff_indices.size());

        //std::cout <<  dof_cell_support[dof].size() << ", " << dof_indices.size() << ", " << dofmatrix_flat.size() << ", " << std::endl;
        
        // Integrate over dofs cell support
        for(std::vector<Cell>::iterator cell = dof_cell_support[dof].begin(); cell != dof_cell_support[dof].end(); ++cell) {
            
            // Extract cell dofs
            for (std::size_t j = 0; j < form_rank; ++j)
              dofs[j] = dofmaps[j]->cell_dofs(cell->index());

            // Create vector of dof-local dofs to cell dofs
            local_dofs.resize(dofs[0].size());
            for (std::size_t j=0; j<dofs[0].size(); j++)
                local_dofs[j] = dofmatrix_to_local_map[dofs[0][j]];
            
            // Find index of current dof in cell dofs
            for (std::size_t q=0; q<dofs[0].size(); q++)
            {
                if (dofs[0][q] == dof)
                {
                    cell_dof_index=q;
                    break;
                }
            }

            // Integrate current cell for each coefficient (with other coeffs set to 0)
            for (coeff_index=0; coeff_index < coeff_indices.size(); coeff_index++)
            {
                // Set 
                _w[coeff_index][cell_dof_index] = 1.0;

                // Update to current cell
                // (UPDATE: Surprisingly expensive, and only vertex coordinates needed. Commented out.)
                //cell->get_cell_data(ufc_cell);
                cell->get_vertex_coordinates(vertex_coordinates);
                //ufc.update(*cell, vertex_coordinates, ufc_cell,
                //        integral->enabled_coefficients());
                
                // Tabulate cell tensor
                //integral->tabulate_tensor(ufc.A.data(), &w_pointer[0],
                //                          vertex_coordinates.data(),
                //                          ufc_cell.orientation);
                integral->tabulate_tensor(ufc.A.data(), &w_pointer[0],
                                          vertex_coordinates.data(),
                                          0);

                for (std::size_t k=0; k<ufc.A.size(); k++)
                {
                    // Fetch dof-local row and column
                    std::size_t row = local_dofs[k/local_dofs.size()];
                    std::size_t column = local_dofs[k%local_dofs.size()];
                    
                    // Get value index in flat dofmatrix
                    dof_flat_index = row*(coeff_indices.size()*dof_indices.size())
                                     +column+coeff_index*dof_indices.size();
                    
                    // Insert into flat dofmatrix
                    dofmatrix_flat[dof_flat_index] += ufc.A[k];

                }
                _w[coeff_index][cell_dof_index] = 0.0;
            }
        }
        
        // Allocate spaces in indices[i] (reclaimed after loop)
        indices[dof].reserve(dof_indices.size()*dof_indices.size());

        // Loop over all rows touched by current dof
        for (std::size_t dof_row=0; dof_row<dof_indices.size(); dof_row++)
        {
            values.clear();
            columns.clear();
            columns.reserve(dof_indices.size());
            values.resize(coeff_indices.size());
            for(coeff_index=0; coeff_index<coeff_indices.size(); coeff_index++)
                values[coeff_index].reserve(dof_indices.size());
            
            for (std::size_t dof_column=0; dof_column<dof_indices.size(); dof_column++)
            {
                for (coeff_index=0; coeff_index<coeff_indices.size(); coeff_index++)
                {
                    dof_flat_index = dof_row*(coeff_indices.size()*dof_indices.size())
                                 +dof_column+coeff_index*dof_indices.size();
                 
                    tmp_values[coeff_index] = dofmatrix_flat[dof_flat_index];   
                }
                
                // Ignore values if all zeros
                if (std::equal(tmp_values.begin(), tmp_values.end(), zero_vector.begin()))
                    continue;
                
                for(coeff_index=0; coeff_index<coeff_indices.size(); coeff_index++)
                    values[coeff_index].push_back(tmp_values[coeff_index]);
                
                columns.push_back(dof_indices[dof_column]);
                
            }

            if (columns.size() == 0)
                continue;
            
            if (columns.size() > _max_cols)
                _max_cols = columns.size();

            // Append index data (row, Ncols, col0, col1, col2, ...)
            row = dof_indices[dof_row];
            indices[dof].push_back(row);
            indices[dof].push_back(columns.size());
            indices[dof].insert(indices[dof].end(), columns.begin(), columns.end());

            // Append values to all_values
            for(coeff_index=0; coeff_index<coeff_indices.size(); coeff_index++)
            {
                all_values.insert(all_values.end(),
                                  values[coeff_index].begin(),
                                  values[coeff_index].end());

            }
        }

        // Reclaim unused memory
        std::vector<dolfin::la_index>(indices[dof]).swap(indices[dof]);

    }

    // Reclaim unused memory
    std::vector<double>(all_values).swap(all_values);
    
    // Print memory usage
    //std::cout << "all_values memory footprint: " << sizeof(double)*all_values.capacity()/(1024*1024) << " MB" << std::endl;
    
    std::size_t indices_capacity=0;
    for (dolfin::la_index i=0; i<N; i++)
    {
        indices_capacity += indices[i].capacity();
    }

    std::cout << "indices memory footprint: " << (sizeof(std::size_t)*indices_capacity)/(1024*1024) << " MB" << std::endl;
    

}


void UnassembledMatrix::assemble(dolfin::GenericTensor& A)
{
    std::vector<std::shared_ptr<const dolfin::Function>> coeffs(_num_coefficients);
    for (std::size_t i=0; i<_num_coefficients; i++)
    {
        coeffs[i] = std::dynamic_pointer_cast<const dolfin::Function>(_a.coefficients()[i]);
        std::cout << (coeffs[i] == NULL) << std::endl;
    }
    assemble(A, coeffs);
}

void UnassembledMatrix::assemble(dolfin::GenericTensor& A, std::vector<std::shared_ptr<const dolfin::Function>> coeffs)
{
    
    if (coeffs.size() != _num_coefficients)
        error("An equal number of coefficients required for assemble as was provided to the constructor.");

    // Init tensor layout
    AssemblerBase().init_global_tensor(A, _a);

    // Down cast to PETScMatrix
    PETScMatrix& _mat = A.down_cast<PETScMatrix>();
    
    //MatSetOption(_mat.mat(), MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
    //MatSetOption(_mat.mat(), MAT_USE_INODES, PETSC_FALSE);
    
    const std::size_t Ndofs=indices.size();
    std::size_t coeff_index;
    
    double * dofvalues = new double[_num_coefficients];
    std::vector<std::vector<double> > local_vector(_num_coefficients);

    // Get local vector of coefficient
    std::size_t offset = dm->ownership_range().first;
    for (coeff_index=0; coeff_index<_num_coefficients; coeff_index++)
    {
        // 1.4.0+:
        std::vector<dolfin::la_index> v(Ndofs);
        std::iota(v.begin(), v.end(), 0);
        local_vector[coeff_index].resize(Ndofs);
        std::cout << "hei-" << coeff_index << std::endl;
        coeffs[coeff_index]->vector()->get_local(local_vector[coeff_index].data(), Ndofs, v.data());
    }
    // Initiate assembled values. No need to resize it in the loop.
    double * assembled_values = new double[_max_cols];
    
    // Set index of all_values-vector
    std::size_t all_values_index = 0;

    // Instantiate indices
    std::size_t idx, dof, indices_idx;
    
    dolfin::la_index row;
    std::size_t Ncols;

    double * coeff_values = new double[_num_coefficients];
    
    // Iterate over all dofs, multiply pre-assembled values with coeff value and add to matrix
    for (std::size_t dof=0; dof<Ndofs; dof++)
    {
        indices_idx = 0;
        while (indices_idx < indices[dof].size())
        {
            row = indices[dof][indices_idx++];
            Ncols = indices[dof][indices_idx++];

            // Assemble values for first coefficient (to avoid resetting to 0.0) before
            // next loop
            for (idx=0; idx<Ncols; idx++)
                assembled_values[idx] = local_vector[0][dof]*all_values[all_values_index++];

            // Loop over all coefficients and assemble into assembled_values
            for (coeff_index=1; coeff_index<_num_coefficients; coeff_index++)
            {
                for (idx=0; idx<Ncols; idx++)
                {
                    assembled_values[idx] += local_vector[coeff_index][dof]*all_values[all_values_index++];
                }
            }
            /*
            std::cout << "Row: " << row << " Ncols: " << Ncols  << std::endl;
            for (std::size_t i=0; i<Ncols; i++){
                std::cout << assembled_values[i] <<", ";
            }
            std::cout << std::endl;
            */
            // Insert single row into matrix
            MatSetValuesLocal(_mat.mat(), 1, &row, Ncols, &indices[dof][indices_idx], assembled_values,
                                          ADD_VALUES);
            

            indices_idx += Ncols;
        }
    }

    // Clean up    
    A.apply("add");
    delete[] assembled_values;
    delete[] dofvalues;
}
