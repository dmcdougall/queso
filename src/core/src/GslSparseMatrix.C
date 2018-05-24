// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


// C++ includes

// Local includes
#include "libmesh/libmesh_config.h"

#include "queso/GslNumericVector.h"
#include "queso/GslSparseMatrix.h"
#include "libmesh/dense_matrix.h"
// #include "libmesh/dof_map.h"
// #include "libmesh/sparsity_pattern.h"

namespace QUESO
{


//-----------------------------------------------------------------------
// GslSparseMatrix members
template <typename T>
void GslSparseMatrix<T>::init (const libMesh::numeric_index_type m_in,
                                 const libMesh::numeric_index_type n_in,
                                 const libMesh::numeric_index_type libmesh_dbg_var(m_l),
                                 const libMesh::numeric_index_type libmesh_dbg_var(n_l),
                                 const libMesh::numeric_index_type nnz,
                                 const libMesh::numeric_index_type,
                                 const libMesh::numeric_index_type)
{
  // noz ignored...  only used for multiple processors!
  libmesh_assert_equal_to (m_in, m_l);
  libmesh_assert_equal_to (n_in, n_l);
  libmesh_assert_equal_to (m_in, n_in);
  libmesh_assert_greater  (nnz, 0);

  this->queso_map.reset(new QUESO::Map(m_in, 0, queso_mpi_comm));
  this->_mat.reset(new QUESO::GslMatrix(queso_env, *queso_map, n_in));

  this->_is_initialized = true;
}



template <typename T>
void GslSparseMatrix<T>::init ()
{
//   // Ignore calls on initialized objects
//   if (this->initialized())
//     return;
//
//   // We need the DofMap for this!
//   libmesh_assert(this->_dof_map);
//
//   // Clear intialized matrices
//   if (this->initialized())
//     this->clear();
//
//   const numeric_index_type n_rows   = this->_dof_map->n_dofs();
//   const numeric_index_type n_cols   = n_rows;
//
// #ifndef NDEBUG
//   // The following variables are only used for assertions,
//   // so avoid declaring them when asserts are inactive.
//   const numeric_index_type n_l = this->_dof_map->n_dofs_on_processor(0);
//   const numeric_index_type m_l = n_l;
// #endif
//
//   // Laspack Matrices only work for uniprocessor cases
//   libmesh_assert_equal_to (n_rows, n_cols);
//   libmesh_assert_equal_to (m_l, n_rows);
//   libmesh_assert_equal_to (n_l, n_cols);
//
//   const std::vector<numeric_index_type> & n_nz = this->_dof_map->get_n_nz();
//
// #ifndef NDEBUG
//   // The following variables are only used for assertions,
//   // so avoid declaring them when asserts are inactive.
//   const std::vector<numeric_index_type> & n_oz = this->_dof_map->get_n_oz();
// #endif
//
//   // Make sure the sparsity pattern isn't empty
//   libmesh_assert_equal_to (n_nz.size(), n_l);
//   libmesh_assert_equal_to (n_oz.size(), n_l);
//
//   if (n_rows==0)
//     {
//       _mat.resize(0,0);
//       return;
//     }
//
//   _mat.resize(n_rows,n_cols);
//   _mat.reserve(n_nz);
//
//   this->_is_initialized = true;
//
//   libmesh_assert_equal_to (n_rows, this->m());
//   libmesh_assert_equal_to (n_cols, this->n());
}



template <typename T>
void GslSparseMatrix<T>::add_matrix(const libMesh::DenseMatrix<T> & dm,
                                      const std::vector<libMesh::numeric_index_type> & rows,
                                      const std::vector<libMesh::numeric_index_type> & cols)

{
  libmesh_assert (this->initialized());
  unsigned int n_rows = libMesh::cast_int<unsigned int>(rows.size());
  unsigned int n_cols = libMesh::cast_int<unsigned int>(cols.size());
  libmesh_assert_equal_to (dm.m(), n_rows);
  libmesh_assert_equal_to (dm.n(), n_cols);


  for (unsigned int i=0; i<n_rows; i++)
    for (unsigned int j=0; j<n_cols; j++)
      this->add(rows[i],cols[j],dm(i,j));
}



template <typename T>
void GslSparseMatrix<T>::get_diagonal (libMesh::NumericVector<T> & dest_in) const
{
  GslNumericVector<T> & dest = libMesh::cast_ref<GslNumericVector<T> &>(dest_in);

  for (unsigned int i = 0; i < dest.size(); i++) {
    dest.set(i, (*_mat)(i,i));
  }
}



template <typename T>
void GslSparseMatrix<T>::get_transpose (libMesh::SparseMatrix<T> & dest_in) const
{
  GslSparseMatrix<T> & dest = libMesh::cast_ref<GslSparseMatrix<T> &>(dest_in);

  *(dest._mat) = _mat->transpose();
}



template <typename T>
GslSparseMatrix<T>::GslSparseMatrix (const libMesh::Parallel::Communicator & comm_in) :
  libMesh::SparseMatrix<T>(comm_in),
  queso_env(),
  queso_mpi_comm(queso_env, comm_in.get()),
  _closed (false)
{
}



template <typename T>
GslSparseMatrix<T>::~GslSparseMatrix ()
{
  this->clear ();
}



template <typename T>
void GslSparseMatrix<T>::clear ()
{
  unsigned int num_cols = 0;
  this->queso_map.reset(new Map(0, 0, this->queso_mpi_comm));
  this->_mat.reset(new GslMatrix(this->queso_env, *(this->queso_map), num_cols));

  _closed = false;
  this->_is_initialized = false;
}



template <typename T>
void GslSparseMatrix<T>::zero ()
{
  _mat->cwSet(0.0);
}



template <typename T>
libMesh::numeric_index_type GslSparseMatrix<T>::m () const
{
  libmesh_assert (this->initialized());

  return _mat->numRowsGlobal();
}



template <typename T>
libMesh::numeric_index_type GslSparseMatrix<T>::n () const
{
  libmesh_assert (this->initialized());

  return _mat->numCols();
}



template <typename T>
libMesh::numeric_index_type GslSparseMatrix<T>::row_start () const
{
  return 0;
}



template <typename T>
libMesh::numeric_index_type GslSparseMatrix<T>::row_stop () const
{
  // return this->m();
}



template <typename T>
void GslSparseMatrix<T>::set (const libMesh::numeric_index_type i,
                                const libMesh::numeric_index_type j,
                                const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());

  // _mat.coeffRef(i,j) = value;
}



template <typename T>
void GslSparseMatrix<T>::add (const libMesh::numeric_index_type i,
                                const libMesh::numeric_index_type j,
                                const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());

  // _mat.coeffRef(i,j) += value;
}



template <typename T>
void GslSparseMatrix<T>::add_matrix(const libMesh::DenseMatrix<T> & dm,
                                      const std::vector<libMesh::numeric_index_type> & dof_indices)
{
  // this->add_matrix (dm, dof_indices, dof_indices);
}



template <typename T>
void GslSparseMatrix<T>::add (const T a_in, libMesh::SparseMatrix<T> & X_in)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (this->m(), X_in.m());
  libmesh_assert_equal_to (this->n(), X_in.n());

  // GslSparseMatrix<T> & X = cast_ref<GslSparseMatrix<T> &> (X_in);

  // _mat += X._mat*a_in;
}




template <typename T>
T GslSparseMatrix<T>::operator () (const libMesh::numeric_index_type i,
                                     const libMesh::numeric_index_type j) const
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());

  // return _mat.coeff(i,j);
}



template <typename T>
libMesh::Real GslSparseMatrix<T>::l1_norm () const
{
  // There does not seem to be a straightforward way to iterate over
  // the columns of an GslSparseMatrix.  So we use some extra
  // storage and keep track of the column sums while going over the
  // row entries...
  std::vector<libMesh::Real> abs_col_sums(this->n());

  // For a row-major Eigen SparseMatrix like we're using, the
  // InnerIterator iterates over the non-zero entries of rows.
  // for (unsigned row=0; row<this->m(); ++row)
  //   {
  //     EigenSM::InnerIterator it(_mat, row);
  //     for (; it; ++it)
  //       abs_col_sums[it.col()] += std::abs(it.value());
  //   }

  return *(std::max_element(abs_col_sums.begin(), abs_col_sums.end()));
}



template <typename T>
libMesh::Real GslSparseMatrix<T>::linfty_norm () const
{
  libMesh::Real max_abs_row_sum = 0.;

  // // For a row-major Eigen SparseMatrix like we're using, the
  // // InnerIterator iterates over the non-zero entries of rows.
  // for (unsigned row=0; row<this->m(); ++row)
  //   {
  //     Real current_abs_row_sum = 0.;
  //     EigenSM::InnerIterator it(_mat, row);
  //     for (; it; ++it)
  //       current_abs_row_sum += std::abs(it.value());
  //
  //     max_abs_row_sum = std::max(max_abs_row_sum, current_abs_row_sum);
  //   }

  return max_abs_row_sum;
}


//------------------------------------------------------------------
// Explicit instantiations
template class GslSparseMatrix<libMesh::Number>;

} // namespace QUESO
