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

#ifdef LIBMESH_HAVE_EIGEN

#include "libmesh/eigen_sparse_vector.h"
#include "libmesh/eigen_sparse_matrix.h"
#include "libmesh/dense_matrix.h"
// #include "libmesh/dof_map.h"
// #include "libmesh/sparsity_pattern.h"
#include <Eigen/SparseLU>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>

// Protecc the fingies
using namespace QUESO;

namespace libMesh
{

template <typename T>
EigenSparseMatrix<T>::EigenSparseMatrix(const EigenSparseVector<T> & v) :
  libMesh::SparseMatrix<T>(
      EigenSparseVector<T>::comm_map.emplace(std::make_pair(
          &(v.queso_map->Comm()),
          libMesh::Parallel::Communicator(
            v.queso_map->Comm().Comm()))).first->second),
  queso_env(new EmptyEnvironment()),
  queso_mpi_comm(v.queso_map->Comm()),
  _closed (false)
{
  this->queso_map.reset(new Map(*v.queso_map));
  this->_is_initialized = true;
}

template <typename T>
EigenSparseMatrix<T>::EigenSparseMatrix(const EigenSparseVector<T> & v, double diagValue) :
  libMesh::SparseMatrix<T>(
      EigenSparseVector<T>::comm_map.emplace(std::make_pair(
          &(v.queso_map->Comm()),
          libMesh::Parallel::Communicator(
            v.queso_map->Comm().Comm()))).first->second),
  queso_env(new EmptyEnvironment()),
  queso_mpi_comm(v.queso_map->Comm()),
  _closed (false)
{
  this->queso_map.reset(new Map(*v.queso_map));

  // Hmmm, is the internal matrix the right size?!
  for (unsigned int i = 0; i < v.sizeLocal(); i++) {
    _mat.coeffRef(i,i) = diagValue;
  }

  this->_is_initialized = true;
}

template <typename T>
EigenSparseMatrix<T>::EigenSparseMatrix(const BaseEnvironment & env,
                                        const Map & map,
                                        unsigned int numCols) :
  libMesh::SparseMatrix<T>(
      EigenSparseVector<T>::comm_map.emplace(std::make_pair(
          &(map.Comm()),
          libMesh::Parallel::Communicator(map.Comm().Comm()))).first->second),
  queso_env(new EmptyEnvironment()),
  queso_mpi_comm(map.Comm()),
  _closed (false)
{
  this->queso_map.reset(new Map(map));
  this->_is_initialized = true;
}


template <typename T>
EigenSparseMatrix<T>::EigenSparseMatrix(const EigenSparseMatrix<T> & B) :
  libMesh::SparseMatrix<T>(
      EigenSparseVector<T>::comm_map.emplace(std::make_pair(
          &(B.queso_map->Comm()),
          libMesh::Parallel::Communicator(
            B.queso_map->Comm().Comm()))).first->second),
  queso_env(new EmptyEnvironment()),
  queso_mpi_comm(B.queso_map->Comm()),
  _closed (false)
{
  this->queso_map.reset(new Map(*B.queso_map));
  this->_is_initialized = true;
}

template <typename T>
unsigned int
EigenSparseMatrix<T>::numCols() const
{
  return this->n();
}

template <typename T>
void
EigenSparseMatrix<T>::zeroLower(bool includeDiagonal)
{
  for (unsigned int i = 1; i < _mat.rows(); i++) {
    for (unsigned int j = 0; j < _mat.cols(); j++) {
      if (j < i) {
        this->set(i,j,0.0);
      }
    }
  }

  if (includeDiagonal) {
    unsigned int num_elements = std::min(_mat.rows(), _mat.cols());
    for (unsigned int i = 0; i < num_elements; i++) {
      this->set(i,i,0.0);
    }
  }
}

template <typename T>
void
EigenSparseMatrix<T>::zeroUpper(bool includeDiagonal)
{
  for (unsigned int i = 0; i < _mat.rows(); i++) {
    for (unsigned int j = 1; j < _mat.cols(); j++) {
      if (j > i) {
        this->set(i,j,0.0);
      }
    }
  }

  if (includeDiagonal) {
    unsigned int num_elements = std::min(_mat.rows(), _mat.cols());
    for (unsigned int i = 0; i < num_elements; i++) {
      this->set(i,i,0.0);
    }
  }
}

template <typename T>
double &
EigenSparseMatrix<T>::operator()(unsigned int i, unsigned int j)
{
  return _mat.coeffRef(i,j);
}

template <typename T>
void
EigenSparseMatrix<T>::cwSet(unsigned int initialTargetRowId,
                            unsigned int initialTargetColId,
                            const EigenSparseMatrix<T> & mat)
{
  queso_require_less_msg(initialTargetRowId, this->numRowsLocal(), "invalid initialTargetRowId");
  queso_require_less_equal_msg((initialTargetRowId + mat.numRowsLocal()), this->numRowsLocal(), "invalid vec.numRowsLocal()");
  queso_require_less_msg(initialTargetColId, this->numCols(), "invalid initialTargetColId");
  queso_require_less_equal_msg((initialTargetColId + mat.numCols()), this->numCols(), "invalid vec.numCols()");

  for (unsigned int i = 0; i < mat.numRowsLocal(); ++i) {
    for (unsigned int j = 0; j < mat.numCols(); ++j) {
      (*this)(initialTargetRowId+i,initialTargetColId+j) = mat(i,j);
    }
  }
}

template <typename T>
unsigned int
EigenSparseMatrix<T>::numRowsLocal() const
{
  return this->n();
}

template <typename T>
EigenSparseMatrix<T> &
EigenSparseMatrix<T>::operator+=(const EigenSparseMatrix<T> & rhs)
{
  add(1.0, const_cast<EigenSparseMatrix<T> &>(rhs));
}

template <typename T>
void
EigenSparseMatrix<T>::mpiSum( const MpiComm & comm, EigenSparseMatrix<T> & M_global) const
{
  // Sanity Checks
  queso_require_equal_to_msg(this->numRowsLocal(), M_global.numRowsLocal(), "local matrices incompatible");
  queso_require_equal_to_msg(this->numCols(), M_global.numCols(), "global matrices incompatible");

  /* TODO: Probably a better way to handle this unpacking/packing of data */
  int size = M_global.numRowsLocal()*M_global.numCols();
  std::vector<double> local( size, 0.0 );
  std::vector<double> global( size, 0.0 );

  int k;
  for( unsigned int i = 0; i < this->numRowsLocal(); i++ ) {
    for( unsigned int j = 0; j < this->numCols(); j++ ) {
      k = i + j*M_global.numCols();
      local[k] = (*this)(i,j);
    }
  }

  comm.Allreduce<double>(&local[0], &global[0], size, RawValue_MPI_SUM,
                 "GslMatrix::mpiSum()",
                 "failed MPI.Allreduce()");

  for( unsigned int i = 0; i < this->numRowsLocal(); i++ ) {
    for( unsigned int j = 0; j < this->numCols(); j++ ) {
      k = i + j*M_global.numCols();
      M_global(i,j) = global[k];
    }
  }
}

template <typename T>
EigenSparseMatrix<T> &
EigenSparseMatrix<T>::operator=(const EigenSparseMatrix<T> & rhs)
{
  queso_mpi_comm = rhs.queso_map->Comm();
  this->queso_map.reset(new Map(*rhs.queso_map));
  _mat = rhs._mat;
  return *this;
}

template <typename T>
void
EigenSparseMatrix<T>::invertMultiply(const EigenSparseMatrix<T> & B, EigenSparseMatrix<T> & X) const
{
  queso_require_equal_to_msg(B.numRowsLocal(), X.numRowsLocal(),
                 "Matrices B and X are incompatible");
  queso_require_equal_to_msg(B.numCols(),      X.numCols(),
                 "Matrices B and X are incompatible");
  queso_require_equal_to_msg(this->numRowsLocal(), X.numRowsLocal(),
                             "This and X matrices are incompatible");

  // Expensive copies -- really bad
  // Need them because col major is required by the solvers in eigen
  Eigen::SparseMatrix<double, Eigen::ColMajor, libMesh::eigen_idx_type> tmpmat = _mat;
  Eigen::SparseMatrix<double, Eigen::ColMajor, libMesh::eigen_idx_type> tmprhs = B._mat;
  Eigen::SparseMatrix<double, Eigen::ColMajor, libMesh::eigen_idx_type> tmpsol = X._mat;

  tmpmat.makeCompressed();
  Eigen::SparseLU<libMesh::EigenSM> solver;
  solver.analyzePattern(tmpmat);
  solver.factorize(tmpmat);
  queso_require_equal_to_msg(solver.info(), Eigen::Success, "decomp failed");

  tmpsol = solver.solve(tmprhs);
  queso_require_equal_to_msg(solver.info(), Eigen::Success, "solve failed");

  // Copy back
  X._mat = tmpsol;
}

template <typename T>
void
EigenSparseMatrix<T>::largestEigen(double & eigenValue, EigenSparseVector<T> & eigenVector) const
{
  unsigned int n = eigenVector.sizeLocal();
  queso_require_not_equal_to_msg(n, 0, "invalid input vector size");

  /*
   * The following notation is used:
   * z = vector used in iteration that ends up being the eigenvector
   *     corresponding to the largest eigenvalue
   * w = vector used in iteration that we extract the largest eigenvalue from.
   */

  // Some parameters associated with the algorithm
  // TODO: Do we want to add the ability to have these set by the user?
  const unsigned int max_num_iterations = 10000;
  const double tolerance = 1.0e-13;

  // Create temporary working vectors.
  // TODO: We shouldn't have to use these - we ought to be able to work directly
  // TODO: with eigenValue and eigenVector. Optimize this once we have regression
  // TODO: tests.
  EigenSparseVector<T> z(*queso_env, *queso_map); // Needs to be initialized to 1.0
  z.cwSet(1.0);
  EigenSparseVector<T> w(*queso_env, *queso_map);

  // Some variables we use inside the loop.
  int index;
  double residual;
  double lambda;

  for( unsigned int k = 0; k < max_num_iterations; ++k) {
    w = (*this) * z;

    // For this algorithm, it's crucial to get the maximum in
    // absolute value, but then to normalize by the actual value
    // and *not* the absolute value.
    EigenSparseVector<T> w_abs(w);
    w_abs.abs();
    index = w_abs.getMaxValueIndex();

    lambda = w[index];

    z = (1.0/lambda) * w;

    // Here we use the norm of the residual as our convergence check:
    // norm( A*x - \lambda*x )
    residual = ( (*this)*z - lambda*z ).norm2();

    if( residual < tolerance ) {
      eigenValue = lambda;

      // TODO: Do we want to normalize this so eigenVector.norm2() = 1?
      eigenVector = z;
      return;
    }
  }

  // If we reach this point, then we didn't converge. Print error message
  // to this effect.
  // TODO: We know we failed at this point - need a way to not need a test
  // TODO: Just specify the error.
  queso_require_less_msg(residual, tolerance, "Maximum num iterations exceeded");
}

template <typename T>
unsigned int
EigenSparseMatrix<T>::numRowsGlobal() const
{
  return m();
}

template <typename T>
double
EigenSparseMatrix<T>::lnDeterminant() const
{
  // Expensive copies -- really bad
  // Need them because col major is required by the solvers in eigen
  Eigen::SparseMatrix<double, Eigen::ColMajor, libMesh::eigen_idx_type> tmpmat = _mat;

  tmpmat.makeCompressed();
  Eigen::SparseLU<libMesh::EigenSM> solver;
  solver.analyzePattern(tmpmat);
  solver.factorize(tmpmat);
  queso_require_equal_to_msg(solver.info(), Eigen::Success, "decomp failed");

  double lndet = solver.logAbsDeterminant();
  queso_require_equal_to_msg(solver.info(), Eigen::Success, "solve failed");

  return lndet;
}

template <typename T>
EigenSparseVector<T>
EigenSparseMatrix<T>::invertMultiply(const EigenSparseVector<T> & b) const
{
  queso_require_equal_to_msg(this->numCols(), b.sizeLocal(), "matrix and rhs have incompatible sizes");

  // Expensive copies -- really bad
  // Need them because col major is required by the solvers in eigen
  Eigen::SparseMatrix<double, Eigen::ColMajor, libMesh::eigen_idx_type> tmpmat = _mat;
  Eigen::Matrix<double, Eigen::Dynamic, 1> tmprhs = b.vec();
  Eigen::Matrix<double, Eigen::Dynamic, 1> tmpsol(m());

  tmpmat.makeCompressed();
  Eigen::SparseLU<libMesh::EigenSM> solver;
  solver.analyzePattern(tmpmat);
  solver.factorize(tmpmat);
  queso_require_equal_to_msg(solver.info(), Eigen::Success, "decomp failed");

  tmpsol = solver.solve(tmprhs);
  queso_require_equal_to_msg(solver.info(), Eigen::Success, "solve failed");

  // Copy back
  EigenSparseVector<T> sol(this->comm(), m(), m());
  sol.vec() = tmpsol;
}

template <typename T>
void
EigenSparseMatrix<T>::invertMultiply(const EigenSparseVector<T> & b, EigenSparseVector<T> & x) const
{
  queso_require_equal_to_msg(this->numCols(), b.sizeLocal(), "matrix and rhs have incompatible sizes");

  // Expensive copies -- really bad
  // Need them because col major is required by the solvers in eigen
  Eigen::SparseMatrix<double, Eigen::ColMajor, libMesh::eigen_idx_type> tmpmat = _mat;
  Eigen::Matrix<double, Eigen::Dynamic, 1> tmprhs = b.vec();
  Eigen::Matrix<double, Eigen::Dynamic, 1> tmpsol = x.vec();

  tmpmat.makeCompressed();
  Eigen::SparseLU<libMesh::EigenSM> solver;
  solver.analyzePattern(tmpmat);
  solver.factorize(tmpmat);
  queso_require_equal_to_msg(solver.info(), Eigen::Success, "decomp failed");

  tmpsol = solver.solve(tmprhs);
  queso_require_equal_to_msg(solver.info(), Eigen::Success, "solve failed");

  // Copy back
  x.vec() = tmpsol;
}


template <typename T>
EigenSparseMatrix<T> &
EigenSparseMatrix<T>::operator*=(double a)
{
  _mat *= a;
  return *this;
}

template <typename T>
int
EigenSparseMatrix<T>::chol()
{
  // Expensive copies -- really bad
  // Need them because col major is required by the solvers in eigen
  Eigen::SparseMatrix<double, Eigen::ColMajor, libMesh::eigen_idx_type> tmpmat = _mat;

  tmpmat.makeCompressed();
  Eigen::SimplicialLLT<libMesh::EigenSM> solver;
  solver.analyzePattern(tmpmat);
  solver.factorize(tmpmat);
  
  if (solver.info() != Eigen::Success) {
    return 1;
  }
  else {
    return 0;
  }
}

template <typename T>
int
EigenSparseMatrix<T>::svd(EigenSparseMatrix<T> & matU,
                          EigenSparseVector<T> & vecS,
                          EigenSparseMatrix & matVt) const
{
  // Expensive copies -- really bad
  // Need them because col major is required by the solvers in eigen
  Eigen::MatrixXd m = _mat;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
  vecS.vec() = svd.singularValues();

  // Sometimes I wonder whether what I'm doing is really the right appraoch.
  // It seems like we should either really be dealing with a dense matrix
  // internally, or abandoning svd and chol and replacing them with iterative
  // linear solvers that are allowed to fail.
  matU.mat() = svd.matrixU().sparseView();
  matVt.mat() = svd.matrixV().sparseView().transpose().eval();
}

template <typename T>
EigenSparseMatrix<T>
EigenSparseMatrix<T>::transpose() const
{
  EigenSparseMatrix<T> answer(*this);
  this->get_transpose(answer);
  return answer;
}

template <typename T>
double
EigenSparseMatrix<T>::determinant() const
{
  return std::exp(this->lnDeterminant());
}

template <typename T>
unsigned int
EigenSparseMatrix<T>::rank(double absoluteZeroThreshold, double relativeZeroThreshold) const
{
  // Expensive copies -- really bad
  // Need them because col major is required by the solvers in eigen
  Eigen::MatrixXd m = _mat;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd relative_values = svd.singularValues();

  if (relative_values[0] > 0.0) {
    relative_values = (1.0 / relative_values[0]) * relative_values;
  }

  unsigned int rank = 0;
  for (unsigned int i = 0; i < relative_values.size(); i++) {
    if ((svd.singularValues()[i] >= absoluteZeroThreshold) &&
        (relative_values[i] >= relativeZeroThreshold)) {
      rank++;
    }
  }

  return rank;
}

template <typename T>
EigenSparseMatrix<T>
EigenSparseMatrix<T>::inverse() const
{
  unsigned int nRows = m();
  unsigned int nCols = n();
  queso_require_equal_to_msg(nRows, nCols, "matrix is not square");

  // Expensive copies -- really bad
  // Need them because col major is required by the solvers in eigen
  Eigen::SparseMatrix<double, Eigen::ColMajor, libMesh::eigen_idx_type> tmpmat = _mat;
  Eigen::SparseMatrix<double, Eigen::ColMajor, libMesh::eigen_idx_type> tmprhs = _mat;
  Eigen::SparseMatrix<double, Eigen::ColMajor, libMesh::eigen_idx_type> tmpsol = _mat;
  tmprhs.setIdentity();

  tmpmat.makeCompressed();
  Eigen::SparseLU<libMesh::EigenSM> solver;
  solver.analyzePattern(tmpmat);
  solver.factorize(tmpmat);
  queso_require_equal_to_msg(solver.info(), Eigen::Success, "decomp failed");

  tmpsol = solver.solve(tmprhs);
  queso_require_equal_to_msg(solver.info(), Eigen::Success, "solve failed");

  // Copy back
  EigenSparseMatrix<T> answer(*this);
  answer._mat = tmpsol;
}

template <typename T>
EigenSparseVector<T>
EigenSparseMatrix<T>::getColumn(const unsigned int column_num) const
{
  EigenSparseVector<T> col(this->comm(), this->m());
  for (unsigned int i = 0; i < col.size(); i++) {
    col[i] = _mat.coeff(i, column_num);
  }

  return col;
}

template <typename T>
void
EigenSparseMatrix<T>::eigen(EigenSparseVector<T> & eigenValues,
                            EigenSparseMatrix<T> * eigenVectors) const
{
  int options;
  if (eigenVectors == NULL) {
    options = Eigen::EigenvaluesOnly;
  }
  else {
    options = Eigen::ComputeEigenvectors;
  }

  // Expensive copies -- really bad
  // Need them because col major (and dense) is required by the solvers in eigen
  Eigen::MatrixXd m = _mat;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(m, options);
  
  eigenValues.vec() = solver.eigenvalues();

  // Sometimes I wonder whether what I'm doing is really the right appraoch.
  // It seems like we should either really be dealing with a dense matrix
  // internally, or abandoning svd and chol and replacing them with iterative
  // linear solvers that are allowed to fail.
  if (options == Eigen::ComputeEigenvectors) {
    eigenVectors->mat() = solver.eigenvectors().sparseView();
  }
}

//-----------------------------------------------------------------------
// EigenSparseMatrix members
template <typename T>
void EigenSparseMatrix<T>::init (const numeric_index_type m_in,
                                 const numeric_index_type n_in,
                                 const numeric_index_type libmesh_dbg_var(m_l),
                                 const numeric_index_type libmesh_dbg_var(n_l),
                                 const numeric_index_type nnz,
                                 const numeric_index_type,
                                 const numeric_index_type)
{
  // noz ignored...  only used for multiple processors!
  libmesh_assert_equal_to (m_in, m_l);
  libmesh_assert_equal_to (n_in, n_l);
  libmesh_assert_equal_to (m_in, n_in);
  libmesh_assert_greater  (nnz, 0);

  _mat.resize(m_in, n_in);
  _mat.reserve(Eigen::Matrix<numeric_index_type, Eigen::Dynamic, 1>::Constant(m_in,nnz));

  this->_is_initialized = true;
}



template <typename T>
void EigenSparseMatrix<T>::init ()
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
void EigenSparseMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
                                      const std::vector<numeric_index_type> & rows,
                                      const std::vector<numeric_index_type> & cols)

{
  libmesh_assert (this->initialized());
  unsigned int n_rows = cast_int<unsigned int>(rows.size());
  unsigned int n_cols = cast_int<unsigned int>(cols.size());
  libmesh_assert_equal_to (dm.m(), n_rows);
  libmesh_assert_equal_to (dm.n(), n_cols);


  for (unsigned int i=0; i<n_rows; i++)
    for (unsigned int j=0; j<n_cols; j++)
      this->add(rows[i],cols[j],dm(i,j));
}



template <typename T>
void EigenSparseMatrix<T>::get_diagonal (NumericVector<T> & dest_in) const
{
  EigenSparseVector<T> & dest = cast_ref<EigenSparseVector<T> &>(dest_in);

  dest._vec = _mat.diagonal();
}



template <typename T>
void EigenSparseMatrix<T>::get_transpose (SparseMatrix<T> & dest_in) const
{
  EigenSparseMatrix<T> & dest = cast_ref<EigenSparseMatrix<T> &>(dest_in);

  dest._mat = _mat.transpose();
}



template <typename T>
EigenSparseMatrix<T>::EigenSparseMatrix (const Parallel::Communicator & comm_in) :
  SparseMatrix<T>(comm_in),
  queso_env(new EmptyEnvironment()),
  queso_mpi_comm(*queso_env, comm_in.get()),
  _closed (false)
{
}



template <typename T>
EigenSparseMatrix<T>::~EigenSparseMatrix ()
{
  this->clear ();
}



template <typename T>
void EigenSparseMatrix<T>::clear ()
{
  _mat.resize(0,0);

  _closed = false;
  this->_is_initialized = false;
}



template <typename T>
void EigenSparseMatrix<T>::zero ()
{
  _mat.setZero();
}



template <typename T>
numeric_index_type EigenSparseMatrix<T>::m () const
{
  libmesh_assert (this->initialized());

  return _mat.rows();
}



template <typename T>
numeric_index_type EigenSparseMatrix<T>::n () const
{
  libmesh_assert (this->initialized());

  return _mat.cols();
}



template <typename T>
numeric_index_type EigenSparseMatrix<T>::row_start () const
{
  return 0;
}



template <typename T>
numeric_index_type EigenSparseMatrix<T>::row_stop () const
{
  return this->m();
}



template <typename T>
void EigenSparseMatrix<T>::set (const numeric_index_type i,
                                const numeric_index_type j,
                                const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());

  _mat.coeffRef(i,j) = value;
}



template <typename T>
void EigenSparseMatrix<T>::add (const numeric_index_type i,
                                const numeric_index_type j,
                                const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());

  _mat.coeffRef(i,j) += value;
}



template <typename T>
void EigenSparseMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
                                      const std::vector<numeric_index_type> & dof_indices)
{
  this->add_matrix (dm, dof_indices, dof_indices);
}



template <typename T>
void EigenSparseMatrix<T>::add (const T a_in, SparseMatrix<T> & X_in)
{
  libmesh_assert (this->initialized());
  libmesh_assert_equal_to (this->m(), X_in.m());
  libmesh_assert_equal_to (this->n(), X_in.n());

  EigenSparseMatrix<T> & X = cast_ref<EigenSparseMatrix<T> &> (X_in);

  _mat += X._mat*a_in;
}




template <typename T>
T EigenSparseMatrix<T>::operator () (const numeric_index_type i,
                                     const numeric_index_type j) const
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->m());
  libmesh_assert_less (j, this->n());

  return _mat.coeff(i,j);
}



template <typename T>
Real EigenSparseMatrix<T>::l1_norm () const
{
  // There does not seem to be a straightforward way to iterate over
  // the columns of an EigenSparseMatrix.  So we use some extra
  // storage and keep track of the column sums while going over the
  // row entries...
  std::vector<Real> abs_col_sums(this->n());

  // For a row-major Eigen SparseMatrix like we're using, the
  // InnerIterator iterates over the non-zero entries of rows.
  for (unsigned row=0; row<this->m(); ++row)
    {
      EigenSM::InnerIterator it(_mat, row);
      for (; it; ++it)
        abs_col_sums[it.col()] += std::abs(it.value());
    }

  return *(std::max_element(abs_col_sums.begin(), abs_col_sums.end()));
}



template <typename T>
Real EigenSparseMatrix<T>::linfty_norm () const
{
  Real max_abs_row_sum = 0.;

  // For a row-major Eigen SparseMatrix like we're using, the
  // InnerIterator iterates over the non-zero entries of rows.
  for (unsigned row=0; row<this->m(); ++row)
    {
      Real current_abs_row_sum = 0.;
      EigenSM::InnerIterator it(_mat, row);
      for (; it; ++it)
        current_abs_row_sum += std::abs(it.value());

      max_abs_row_sum = std::max(max_abs_row_sum, current_abs_row_sum);
    }

  return max_abs_row_sum;
}


//------------------------------------------------------------------
// Explicit instantiations
template class EigenSparseMatrix<Number>;

} // namespace libMesh

namespace QUESO
{

template <typename T>
libMesh::EigenSparseVector<T> operator*(const libMesh::EigenSparseMatrix<T> & mat,
                                        const libMesh::EigenSparseVector<T> & vec)
{
  libMesh::EigenSparseVector<T> answer(vec.comm(), mat.m());
  answer.vec() = (mat.mat() * vec.vec()).eval();
  return answer;
}

template <typename T>
libMesh::EigenSparseMatrix<T> operator*(double a, const libMesh::EigenSparseMatrix<T> & mat)
{
  libMesh::EigenSparseMatrix<T> answer(mat);
  answer.mat() *= a;
  return answer;
}

template <typename T>
libMesh::EigenSparseMatrix<T> operator*(const libMesh::EigenSparseMatrix<T> & m1,
                                        const libMesh::EigenSparseMatrix<T> & m2)
{
  libMesh::EigenSparseMatrix<T> answer(m1.comm());
  answer.init(m1.m(), m2.n(), m1.m(), m2.n());
  answer.mat() = m1.mat() * m2.mat();
  return answer;
}

template <typename T>
libMesh::EigenSparseMatrix<T> matrixProduct(const libMesh::EigenSparseVector<T> & v1, const libMesh::EigenSparseVector<T> & v2)
{
  unsigned int nRows = v1.sizeLocal();
  unsigned int nCols = v2.sizeLocal();
  libMesh::EigenSparseMatrix<T> answer(v1.comm());
  answer.init(nRows, nCols, nRows, nCols);

  for (unsigned int i = 0; i < nRows; ++i) {
    double value1 = v1[i];
    for (unsigned int j = 0; j < nCols; ++j) {
      answer(i,j) = value1*v2[j];
    }
  }

  return answer;
}

template <typename T>
libMesh::EigenSparseMatrix<T> operator+(const libMesh::EigenSparseMatrix<T> & m1,
                                        const libMesh::EigenSparseMatrix<T> & m2)
{
  libMesh::EigenSparseMatrix<T> mat(m1);
  mat.mat() += m2.mat();
  return mat;
}

template libMesh::EigenSparseMatrix<libMesh::Number> operator*(double, const libMesh::EigenSparseMatrix<libMesh::Number> &);
template libMesh::EigenSparseMatrix<libMesh::Number> operator*(const libMesh::EigenSparseMatrix<libMesh::Number> &, const libMesh::EigenSparseMatrix<libMesh::Number> &);
template libMesh::EigenSparseMatrix<libMesh::Number> operator+(const libMesh::EigenSparseMatrix<libMesh::Number> &, const libMesh::EigenSparseMatrix<libMesh::Number> &);
template libMesh::EigenSparseMatrix<libMesh::Number> matrixProduct(const libMesh::EigenSparseVector<libMesh::Number> &, const libMesh::EigenSparseVector<libMesh::Number> &);

}

#endif // #ifdef LIBMESH_HAVE_EIGEN
