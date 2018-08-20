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
#include "queso/GslNumericBlockMatrix.h"
#include "libmesh/dense_matrix.h"
// #include "libmesh/dof_map.h"
// #include "libmesh/sparsity_pattern.h"
#include <queso/GslSparseMatrix.h>

namespace QUESO
{


//-----------------------------------------------------------------------
// GslNumericBlockMatrix members
template <typename T>
void GslNumericBlockMatrix<T>::init (const libMesh::numeric_index_type m_in,
                                 const libMesh::numeric_index_type n_in,
                                 const libMesh::numeric_index_type libmesh_dbg_var(m_l),
                                 const libMesh::numeric_index_type libmesh_dbg_var(n_l),
                                 const libMesh::numeric_index_type nnz,
                                 const libMesh::numeric_index_type,
                                 const libMesh::numeric_index_type)
{
  queso_not_implemented();
  // // noz ignored...  only used for multiple processors!
  // libmesh_assert_equal_to (m_in, m_l);
  // libmesh_assert_equal_to (n_in, n_l);
  // libmesh_assert_equal_to (m_in, n_in);
  // libmesh_assert_greater  (nnz, 0);
  //
  // this->queso_map.reset(new QUESO::Map(m_in, 0, queso_mpi_comm));
  // this->_mat.reset(new QUESO::GslMatrix(*queso_env, *queso_map, n_in));
  //
  // this->_is_initialized = true;
}



template <typename T>
void GslNumericBlockMatrix<T>::init ()
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
void GslNumericBlockMatrix<T>::add_matrix(const libMesh::DenseMatrix<T> & dm,
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
void GslNumericBlockMatrix<T>::get_diagonal (libMesh::NumericVector<T> & dest_in) const
{
  queso_not_implemented();
  // GslNumericVector<T> & dest = libMesh::cast_ref<GslNumericVector<T> &>(dest_in);
  //
  // for (unsigned int i = 0; i < dest.size(); i++) {
  //   dest.set(i, (*_mat)(i,i));
  // }
}



template <typename T>
void GslNumericBlockMatrix<T>::get_transpose (libMesh::SparseMatrix<T> & dest_in) const
{
  queso_not_implemented();
  // GslNumericBlockMatrix<T> & dest = libMesh::cast_ref<GslNumericBlockMatrix<T> &>(dest_in);
  //
  // *(dest._mat) = _mat->transpose();
}



// template <typename T>
// GslNumericBlockMatrix<T>::GslNumericBlockMatrix (const libMesh::Parallel::Communicator & comm_in) :
//   libMesh::SparseMatrix<T>(comm_in),
//   queso_env(new EmptyEnvironment()),
//   queso_mpi_comm(*queso_env, comm_in.get()),
//   _closed (false)
// {
// }



template <typename T>
GslNumericBlockMatrix<T>::~GslNumericBlockMatrix ()
{
  this->clear ();
}

// template <typename T>
// GslNumericBlockMatrix<T>::GslNumericBlockMatrix(const GslNumericVector<T> & v) :
//   libMesh::SparseMatrix<T>(
//       GslNumericVector<T>::comm_map.emplace(std::make_pair(
//           &(v.queso_map->Comm()),
//           libMesh::Parallel::Communicator(
//             v.queso_map->Comm().Comm()))).first->second),
//   queso_env(new EmptyEnvironment()),
//   queso_mpi_comm(v._vec->map().Comm()),
//   _closed (false)
// {
//   queso_not_implemented();
//   // this->queso_map.reset(new Map(v._vec->map()));
//   // this->_mat.reset(new GslMatrix(v._vec->env(), v._vec->map(), v._vec->sizeLocal()));
// }
//
// template <typename T>
// GslNumericBlockMatrix<T>::GslNumericBlockMatrix(const GslNumericVector<T> & v, double diagValue) :
//   libMesh::SparseMatrix<T>(
//       GslNumericVector<T>::comm_map.emplace(std::make_pair(
//           &(v.queso_map->Comm()),
//           libMesh::Parallel::Communicator(
//             v.queso_map->Comm().Comm()))).first->second),
//   queso_env(new EmptyEnvironment()),
//   queso_mpi_comm(v._vec->map().Comm()),
//   _closed (false)
// {
//   queso_not_implemented();
//   // this->queso_map.reset(new Map(v._vec->map()));
//   // this->_mat.reset(new GslMatrix(*v._vec, diagValue));
// }

// template <typename T>
// GslNumericBlockMatrix<T>::GslNumericBlockMatrix(const GslNumericBlockMatrix<T> & B) :
//   libMesh::SparseMatrix<T>(
//       GslNumericVector<T>::comm_map.emplace(std::make_pair(
//           &(B.queso_map->Comm()),
//           libMesh::Parallel::Communicator(
//             B.queso_map->Comm().Comm()))).first->second),
//   queso_env(new EmptyEnvironment()),
//   queso_mpi_comm(B._mat->map().Comm()),
//   _closed (false)
// {
//   queso_not_implemented();
//   // this->queso_map.reset(new Map(B._mat->map()));
//   // this->_mat.reset(new GslMatrix(*B._mat));
// }


// template <typename T>
// GslNumericBlockMatrix<T>::GslNumericBlockMatrix(const BaseEnvironment & env, const Map & map, unsigned int numCols) :
//   libMesh::SparseMatrix<T>(
//       GslNumericVector<T>::comm_map.emplace(std::make_pair(
//           &(map.Comm()),
//           libMesh::Parallel::Communicator(map.Comm().Comm()))).first->second),
//   queso_env(new EmptyEnvironment()),
//   queso_mpi_comm(map.Comm()),
//   _closed (false)
// {
//   queso_not_implemented();
//   // this->queso_map.reset(new Map(map));
//   // this->_mat.reset(new GslMatrix(env, map, numCols));
// }


template <typename T>
void GslNumericBlockMatrix<T>::clear ()
{
  queso_not_implemented();
  // unsigned int num_cols = 0;
  // this->queso_map.reset(new Map(0, 0, this->queso_mpi_comm));
  // this->_mat.reset(new GslMatrix(*(this->queso_env), *(this->queso_map), num_cols));
  //
  // _closed = false;
  // this->_is_initialized = false;
}



template <typename T>
void GslNumericBlockMatrix<T>::zero ()
{
  queso_not_implemented();
  // _mat->cwSet(0.0);
}



template <typename T>
libMesh::numeric_index_type GslNumericBlockMatrix<T>::m () const
{
  libmesh_assert (this->initialized());

  return _mat->numRowsGlobal();
}



template <typename T>
libMesh::numeric_index_type GslNumericBlockMatrix<T>::n () const
{
  libmesh_assert (this->initialized());

  return _mat->numCols();
}



template <typename T>
libMesh::numeric_index_type GslNumericBlockMatrix<T>::row_start () const
{
  return 0;
}



template <typename T>
libMesh::numeric_index_type GslNumericBlockMatrix<T>::row_stop () const
{
  return this->m();
}



template <typename T>
void GslNumericBlockMatrix<T>::set (const libMesh::numeric_index_type i,
                                const libMesh::numeric_index_type j,
                                const T value)
{
  queso_not_implemented();
  // libmesh_assert (this->initialized());
  // libmesh_assert_less (i, this->m());
  // libmesh_assert_less (j, this->n());
  //
  // (*_mat)(i,j) = value;
}



template <typename T>
void GslNumericBlockMatrix<T>::add (const libMesh::numeric_index_type i,
                                const libMesh::numeric_index_type j,
                                const T value)
{
  queso_not_implemented();
  // libmesh_assert (this->initialized());
  // libmesh_assert_less (i, this->m());
  // libmesh_assert_less (j, this->n());
  //
  // (*_mat)(i,j) += value;
}



template <typename T>
void GslNumericBlockMatrix<T>::add_matrix(const libMesh::DenseMatrix<T> & dm,
                                      const std::vector<libMesh::numeric_index_type> & dof_indices)
{
  this->add_matrix (dm, dof_indices, dof_indices);
}



template <typename T>
void GslNumericBlockMatrix<T>::add (const T a_in, libMesh::SparseMatrix<T> & X_in)
{
  queso_not_implemented();
  // libmesh_assert (this->initialized());
  // libmesh_assert_equal_to (this->m(), X_in.m());
  // libmesh_assert_equal_to (this->n(), X_in.n());
  //
  // GslNumericBlockMatrix<T> & X = libMesh::cast_ref<GslNumericBlockMatrix<T> &> (X_in);
  //
  // // We don't guarantee X_in is preserved by this function, so we can modify it
  // // while scaling by a_in here.
  // *(X._mat) *= a_in;
  //
  // *_mat += *(X._mat);
}




template <typename T>
T GslNumericBlockMatrix<T>::operator () (const libMesh::numeric_index_type i,
                                     const libMesh::numeric_index_type j) const
{
  queso_not_implemented();
  // libmesh_assert (this->initialized());
  // libmesh_assert_less (i, this->m());
  // libmesh_assert_less (j, this->n());
  //
  // return (*_mat)(i,j);
}



template <typename T>
libMesh::Real GslNumericBlockMatrix<T>::l1_norm () const
{
  queso_not_implemented();
  // double l1_norm = 0.0;
  //
  // GslVector col_vec(*queso_env, *queso_map);
  //
  // for (unsigned col=0; col<this->n(); ++col) {
  //   _mat->getColumn(col, col_vec);
  //   l1_norm = std::max(l1_norm, col_vec.norm1());
  // }
  //
  // return l1_norm;
}



template <typename T>
libMesh::Real GslNumericBlockMatrix<T>::linfty_norm () const
{
  queso_not_implemented();
  // double linfty_norm = 0.0;
  //
  // GslVector row_vec(*queso_env, *queso_map);
  //
  // for (unsigned row=0; row<this->m(); ++row) {
  //   _mat->getRow(row, row_vec);
  //   linfty_norm = std::max(linfty_norm, row_vec.norm1());
  // }
  //
  // return linfty_norm;
}

template <typename T>
unsigned int
GslNumericBlockMatrix<T>::numCols() const
{
  return this->_mat->numCols();
}

template <typename T>
void
GslNumericBlockMatrix<T>::zeroLower(bool includeDiagonal)
{
  return this->_mat->zeroLower(includeDiagonal);
}

template <typename T>
void
GslNumericBlockMatrix<T>::zeroUpper(bool includeDiagonal)
{
  return this->_mat->zeroUpper(includeDiagonal);
}

// template <typename T>
// double &
// GslNumericBlockMatrix<T>::operator()(unsigned int i, unsigned int j)
// {
//   queso_not_implemented();
//   // return (*_mat)(i, j);
// }

// template <typename T>
// void
// GslNumericBlockMatrix<T>::cwSet(unsigned int rowId, unsigned int colId, const GslNumericBlockMatrix<T> & mat)
// {
//   this->_mat->cwSet(rowId, colId, *mat._mat);
// }

template <typename T>
int
GslNumericBlockMatrix<T>::chol()
{
  return this->_mat->chol();
}

// template <typename T>
// int
// GslNumericBlockMatrix<T>::svd(GslNumericBlockMatrix<T> & matU,
//                         GslNumericVector<T> & vecS,
//                         GslNumericBlockMatrix & matVt) const
// {
//   queso_not_implemented();
//   // return this->_mat->svd(*matU._mat, *vecS._vec, *matVt._mat);
// }

// template <typename T>
// GslNumericVector<T>
// GslNumericBlockMatrix<T>::invertMultiply(const GslNumericVector<T> & b) const
// {
//   GslNumericVector<T> answer(b);
//   *(answer._vec) = this->_mat->invertMultiply(*b._vec);
//   return answer;
// }

// template <typename T>
// void
// GslNumericBlockMatrix<T>::invertMultiply(const GslNumericBlockMatrix<T> & B, GslNumericBlockMatrix<T> X) const
// {
//   this->_mat->invertMultiply(*B._mat, *X._mat);
// }

template <typename T>
void
GslNumericBlockMatrix<T>::invertMultiply(const GslNumericVector<T> & b, GslNumericVector<T> & x) const
{
  this->_mat->invertMultiply(*b._vec, *x._vec);
}

//
// template <typename T>
// double
// GslNumericBlockMatrix<T>::lnDeterminant() const
// {
//   return this->_mat->lnDeterminant();
// }

template <typename T>
unsigned int
GslNumericBlockMatrix<T>::numRowsGlobal() const
{
  return this->_mat->numRowsGlobal();
}

// template <typename T>
// GslNumericBlockMatrix<T> &
// GslNumericBlockMatrix<T>::operator=(const GslNumericBlockMatrix<T> & rhs)
// {
//   queso_mpi_comm = rhs._mat->map().Comm();
//   this->queso_map.reset(new Map(rhs._mat->map()));
//   *(this->_mat) = *rhs._mat;
//   return *this;
// }
//
// template <typename T>
// GslNumericBlockMatrix<T> &
// GslNumericBlockMatrix<T>::operator*=(double a)
// {
//   (*this->_mat) *= a;
//   return *this;
// }
//
// template <typename T>
// GslNumericBlockMatrix<T> &
// GslNumericBlockMatrix<T>::operator+=(const GslNumericBlockMatrix<T> & rhs)
// {
//   (*this->_mat) += *rhs._mat;
//   return *this;
// }
//
// template <typename T>
// GslNumericVector<T>
// GslNumericBlockMatrix<T>::multiply(const GslNumericVector<T> & x) const
// {
//   GslNumericVector<T> answer(x);
//   *answer._vec = this->_mat->multiply(*x._vec);
//   return answer;
// }
//
// template <typename T>
// void
// GslNumericBlockMatrix<T>::largestEigen(double & eigenValue, GslNumericVector<T> & eigenVector) const
// {
//   this->_mat->largestEigen(eigenValue, *eigenVector._vec);
// }
//
// template <typename T>
// void
// GslNumericBlockMatrix<T>::mpiSum(const MpiComm & comm, GslNumericBlockMatrix<T> & M_global) const
// {
//   this->_mat->mpiSum(comm, *M_global._mat);
// }
//
// template <typename T>
// GslNumericBlockMatrix<T>
// GslNumericBlockMatrix<T>::transpose() const
// {
//   queso_not_implemented();
//   // GslMatrix internal_answer(_mat->transpose());
//   //
//   // int num_elements = internal_answer.map().NumGlobalElements();
//   // GslNumericBlockMatrix<T> answer(internal_answer.env(),
//   //                           internal_answer.map(),
//   //                           num_elements);
//   //
//   // *(answer._mat) = internal_answer;
//   //
//   // return answer;
// }
//
// template <typename T>
// const Map &
// GslNumericBlockMatrix<T>::map() const
// {
//   return this->_mat->map();
// }
//
// template <typename T>
// void
// GslNumericBlockMatrix<T>::eigen(GslNumericVector<T> & eigenValues, GslNumericBlockMatrix<T> * eigenVectors) const
// {
//   this->_mat->eigen(*eigenValues._vec, eigenVectors->_mat.get());
// }
//
// template <typename T>
// GslNumericVector<T>
// GslNumericBlockMatrix<T>::getColumn(const unsigned int column_num) const
// {
//   GslVector internal_answer(this->_mat->getColumn(column_num));
//   GslNumericVector<T> answer(internal_answer.env(), internal_answer.map());
//
//   *answer._vec = internal_answer;
//
//   return answer;
// }
//
// template <typename T>
// GslNumericBlockMatrix<T>
// GslNumericBlockMatrix<T>::inverse() const
// {
//   GslMatrix internal_answer(this->_mat->inverse());
//
//   GslNumericBlockMatrix<T> answer(internal_answer.env(),
//                             internal_answer.map(),
//                             internal_answer.map().NumGlobalElements());
//
//   return answer;
// }
//
// template <typename T>
// unsigned int
// GslNumericBlockMatrix<T>::rank(double absoluteZeroThreshold, double relativeZeroThreshold) const
// {
//   return this->_mat->rank(absoluteZeroThreshold, relativeZeroThreshold);
// }
//
// template <typename T>
// double
// GslNumericBlockMatrix<T>::determinant() const
// {
//   return this->_mat->determinant();
// }
//
// template <typename T>
// const BaseEnvironment &
// GslNumericBlockMatrix<T>::env() const
// {
//   queso_not_implemented();
//   // return this->_mat->env();
// }
//
// template <typename T>
// GslNumericVector<T> operator*(const GslNumericBlockMatrix<T> & mat,
//                               const GslNumericVector<T> & vec)
// {
//   queso_not_implemented();
//   // return mat.multiply(vec);
// }
//
// template <typename T>
// GslNumericBlockMatrix<T> operator*(double a,
//                              const GslNumericBlockMatrix<T> & mat)
// {
//   queso_not_implemented();
//   // GslNumericBlockMatrix<T> answer(mat);
//   // answer *= a;
//   // return answer;
// }
//
// template <typename T>
// GslNumericBlockMatrix<T> matrixProduct(const GslNumericVector<T> & v1,
//                                  const GslNumericVector<T> & v2)
// {
//   queso_not_implemented();
//   // unsigned int nRows = v1.sizeLocal();
//   // unsigned int nCols = v2.sizeLocal();
//   // GslNumericBlockMatrix<T> answer(v1.env(),v1.map(),nCols);
//   //
//   // for (unsigned int i = 0; i < nRows; ++i) {
//   //   double value1 = v1[i];
//   //   for (unsigned int j = 0; j < nCols; ++j) {
//   //     answer(i,j) = value1*v2[j];
//   //   }
//   // }
//   //
//   // return answer;
// }
//
// template <typename T>
// GslNumericBlockMatrix<T> operator*(const GslNumericBlockMatrix<T> & m1,
//                              const GslNumericBlockMatrix<T> & m2)
// {
//   unsigned int m1Rows = m1.numRowsGlobal();  // Because they're serial
//   unsigned int m1Cols = m1.numCols();
//   unsigned int m2Rows = m2.numRowsGlobal();  // Because they're serial
//   unsigned int m2Cols = m2.numCols();
//
//   queso_require_equal_to_msg(m1Cols, m2Rows, "different sizes m1Cols and m2Rows");
//
//   GslNumericBlockMatrix<T> mat(m1.env(),m1.map(),m2Cols);
//
//   unsigned int commonSize = m1Cols;
//   for (unsigned int row1 = 0; row1 < m1Rows; ++row1) {
//     for (unsigned int col2 = 0; col2 < m2Cols; ++col2) {
//       double result = 0.;
//       for (unsigned int k = 0; k < commonSize; ++k) {
//         result += m1(row1,k)*m2(k,col2);
//       }
//       mat(row1,col2) = result;
//     }
//   }
//
//   return mat;
// }
//
// template <typename T>
// GslNumericBlockMatrix<T> operator+(const GslNumericBlockMatrix<T> & m1,
//                              const GslNumericBlockMatrix<T> & m2)
// {
//   GslNumericBlockMatrix<T> mat(m1);
//   mat += m2;
//   return mat;
// }

template <typename T>
unsigned int
GslNumericBlockMatrix<T>::numRowsLocal() const
{
  return this->_mat->numRowsLocal();
}

// template <typename T>
// void
// GslNumericBlockMatrix<T>::subReadContents(const std::string & fileName,
//                                     const std::string & fileType,
//                                     const std::set<unsigned int> & allowedSubEnvIds)
// {
//   this->_mat->subReadContents(fileName, fileType, allowedSubEnvIds);
// }
//
// template <typename T>
// void
// GslNumericBlockMatrix<T>::cwSet(double value)
// {
//   this->_mat->cwSet(value);
// }
//
// template <typename T>
// void
// GslNumericBlockMatrix<T>::subWriteContents(const std::string & varNamePrefix,
//                                      const std::string & fileName,
//                                      const std::string & fileType,
//                                      const std::set<unsigned int> & allowedSubEnvIds) const
// {
//   this->_mat->subWriteContents(varNamePrefix, fileName, fileType, allowedSubEnvIds);
// }
//
// template <typename T>
// GslNumericBlockMatrix<T> &
// GslNumericBlockMatrix<T>::operator/=(double a)
// {
//   *this->_mat /= a;
//   return *this;
// }

template <typename T>
unsigned int
GslNumericBlockMatrix<T>::numBlocks() const
{
  return this->_mat->numBlocks();
}

template <typename T>
GslSparseMatrix<T> &
GslNumericBlockMatrix<T>::getBlock(unsigned int i) const
{
  GslMatrix internal_answer(this->_mat->getBlock(i));

  GslSparseMatrix<T> answer(internal_answer.env(),
                            internal_answer.map(),
                            internal_answer.map().NumGlobalElements());

  return answer;
}

//------------------------------------------------------------------
// Explicit instantiations
template class GslNumericBlockMatrix<libMesh::Number>;

// template GslNumericVector<libMesh::Number> operator*(const GslNumericBlockMatrix<libMesh::Number> &, const GslNumericVector<libMesh::Number> &);
// template GslNumericBlockMatrix<libMesh::Number> operator*(double, const GslNumericBlockMatrix<libMesh::Number> &);
// template GslNumericBlockMatrix<libMesh::Number> matrixProduct(const GslNumericVector<libMesh::Number> &, const GslNumericVector<libMesh::Number> &);
// template GslNumericBlockMatrix<libMesh::Number> operator*(const GslNumericBlockMatrix<libMesh::Number> &, const GslNumericBlockMatrix<libMesh::Number> &);
// template GslNumericBlockMatrix<libMesh::Number> operator+(const GslNumericBlockMatrix<libMesh::Number> &, const GslNumericBlockMatrix<libMesh::Number> &);

} // namespace QUESO
