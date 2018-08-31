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

#ifndef QUESO_GSL_NUMERIC_BLOCK_MATRIX_H
#define QUESO_GSL_NUMERIC_BLOCK_MATRIX_H

#include "libmesh/libmesh_config.h"

// Local includes
#include "libmesh/sparse_matrix.h"

#include <queso/GslBlockMatrix.h>

// Eigen includes

// C++ includes
#include <algorithm>
#include <cstddef>

namespace libMesh
{
// Forward declarations
template <typename T> class DenseMatrix;
// template <typename T> class EigenSparseLinearSolver;
}

namespace QUESO
{

// Forward declarations
template <typename T> class GslNumericVector;

/**
 * The GslNumericBlockMatrix class wraps a sparse matrix object from the GSL
 * library.
 *
 * \author Damon McDougall
 * \date 2018
 */
template <typename T>
class GslNumericBlockMatrix libmesh_final : public libMesh::SparseMatrix<T>
{

public:
  /**
   * Constructor; initializes the matrix to
   * be empty, without any structure, i.e.
   * the matrix is not usable at all. This
   * constructor is therefore only useful
   * for matrices which are members of a
   * class. All other matrices should be
   * created at a point in the data flow
   * where all necessary information is
   * available.
   *
   * You have to initialize
   * the matrix before usage with
   * \p init(...).
   */
  // GslNumericBlockMatrix (const libMesh::Parallel::Communicator & comm
  //                    LIBMESH_CAN_DEFAULT_TO_COMMWORLD);

  /**
   * Destructor. Free all memory, but do not
   * release the memory of the sparsity
   * structure.
   */
  ~GslNumericBlockMatrix ();

  // Queso ctors and methods
  // GslNumericBlockMatrix(const GslNumericVector<T> & v);
  // GslNumericBlockMatrix(const GslNumericVector<T> & v, double diagValue);
  // GslNumericBlockMatrix(const GslNumericBlockMatrix<T> & B);
  // GslNumericBlockMatrix(const BaseEnvironment & env, const Map & map, unsigned int numCols);
  unsigned int numCols() const;
  void zeroLower(bool includeDiagonal=false);
  void zeroUpper(bool includeDiagonal=false);
  // double & operator()(unsigned int i, unsigned int j);
  // void cwSet(unsigned int rowId, unsigned int colId, const GslNumericBlockMatrix<T> & mat);
  int chol();
  // int svd(GslNumericBlockMatrix<T> & matU,
  //         GslNumericVector<T> & vecS,
  //         GslNumericBlockMatrix & matVt) const;
  // GslNumericVector<T> invertMultiply(const GslNumericVector<T> & b) const;
  // void invertMultiply(const GslNumericBlockMatrix<T> & B, GslNumericBlockMatrix<T> X) const;
  void invertMultiply(const GslNumericVector<T> & b, GslNumericVector<T> & x) const;
  // double lnDeterminant() const;
  unsigned int numRowsGlobal() const;
  // GslNumericBlockMatrix<T> & operator=(const GslNumericBlockMatrix & rhs);
  // GslNumericBlockMatrix<T> & operator*=(double a);
  // GslNumericBlockMatrix<T> & operator+=(const GslNumericBlockMatrix<T> & rhs);
  // GslNumericVector<T> multiply(const GslNumericVector<T> & x) const;
  // void largestEigen(double & eigenValue, GslNumericVector<T> & eigenVector) const;
  // void mpiSum(const MpiComm & comm, GslNumericBlockMatrix<T> & M_global) const;
  // GslNumericBlockMatrix<T> transpose() const;
  // const Map & map() const;
  // void eigen(GslNumericVector<T> & eigenValues, GslNumericBlockMatrix<T> * eigenVectors) const;
  // GslNumericVector<T> getColumn(const unsigned int column_num) const;
  // GslNumericBlockMatrix<T> inverse() const;
  // unsigned int rank(double absoluteZeroThreshold, double relativeZeroThreshold) const;
  // double determinant() const;
  // const BaseEnvironment & env() const;
  unsigned int numRowsLocal() const;
  // void subReadContents(const std::string & fileName,
  //                      const std::string & fileType,
  //                      const std::set<unsigned int> & allowedSubEnvIds);
  // void cwSet(double value);
  // void subWriteContents(const std::string & varNamePrefix,
  //                       const std::string & fileName,
  //                       const std::string & fileType,
  //                       const std::set<unsigned int> & allowedSubEnvIds) const;
  // GslNumericBlockMatrix<T> & operator/=(double a);
  unsigned int numBlocks() const;
  GslSparseMatrix<T> & getBlock(unsigned int i) const;


  /**
   * Convenient typedefs
   */
  typedef std::unique_ptr<GslBlockMatrix> DataType;

  /**
   * Initialize a QUESO GSL matrix that is of global
   * dimension \f$ m \times  n \f$ with local dimensions
   * \f$ m_l \times n_l \f$.  \p nnz is the number of on-processor
   * nonzeros per row (defaults to 30).
   * \p noz is the number of on-processor
   * nonzeros per row (defaults to 30).
   * Optionally supports a block size, which indicates dense coupled blocks
   * for systems with multiple variables all of the same type.
   */
  virtual void init (const libMesh::numeric_index_type m,
                     const libMesh::numeric_index_type n,
                     const libMesh::numeric_index_type m_l,
                     const libMesh::numeric_index_type n_l,
                     const libMesh::numeric_index_type nnz=30,
                     const libMesh::numeric_index_type noz=10,
                     const libMesh::numeric_index_type blocksize=1) libmesh_override;

  /**
   * Initialize using sparsity structure computed by \p dof_map.
   */
  virtual void init () libmesh_override;

  /**
   * Release all memory and return
   * to a state just like after
   * having called the default
   * constructor.
   */
  virtual void clear () libmesh_override;

  /**
   * Set all entries to 0.
   */
  virtual void zero () libmesh_override;

  /**
   * Close the matrix.  Dummy routine.  After calling
   * this method \p closed() is true and the matrix can
   * be used in computations.
   */
  virtual void close () const libmesh_override { const_cast<GslNumericBlockMatrix<T> *>(this)->_closed = true; }

  /**
   * @returns \p m, the row-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */
  virtual libMesh::numeric_index_type m () const libmesh_override;

  /**
   * @returns \p n, the column-dimension of
   * the matrix where the marix is \f$ M \times N \f$.
   */
  virtual libMesh::numeric_index_type n () const libmesh_override;

  /**
   * return row_start, the index of the first
   * matrix row stored on this processor
   */
  virtual libMesh::numeric_index_type row_start () const libmesh_override;

  /**
   * return row_stop, the index of the last
   * matrix row (+1) stored on this processor
   */
  virtual libMesh::numeric_index_type row_stop () const libmesh_override;

  /**
   * Set the element \p (i,j) to \p value.
   * Throws an error if the entry does
   * not exist. Still, it is allowed to store
   * zero values in non-existent fields.
   */
  virtual void set (const libMesh::numeric_index_type i,
                    const libMesh::numeric_index_type j,
                    const T value) libmesh_override;

  /**
   * Add \p value to the element
   * \p (i,j).  Throws an error if
   * the entry does not
   * exist. Still, it is allowed to
   * store zero values in
   * non-existent fields.
   */
  virtual void add (const libMesh::numeric_index_type i,
                    const libMesh::numeric_index_type j,
                    const T value) libmesh_override;

  /**
   * Add the full matrix to the
   * GSL matrix.  This is useful
   * for adding an element matrix
   * at assembly time
   */
  virtual void add_matrix (const libMesh::DenseMatrix<T> & dm,
                           const std::vector<libMesh::numeric_index_type> & rows,
                           const std::vector<libMesh::numeric_index_type> & cols) libmesh_override;

  /**
   * Same, but assumes the row and column maps are the same.
   * Thus the matrix \p dm must be square.
   */
  virtual void add_matrix (const libMesh::DenseMatrix<T> & dm,
                           const std::vector<libMesh::numeric_index_type> & dof_indices) libmesh_override;

  /**
   * Add a Sparse matrix \p X, scaled with \p a, to \p this,
   * stores the result in \p this: \f$\texttt{this} += a*X \f$.
   * \p LASPACK does not provide a true \p axpy for matrices,
   * so a hand-coded version with hopefully acceptable performance
   * is provided.
   */
  virtual void add (const T a, libMesh::SparseMatrix<T> & X) libmesh_override;

  /**
   * Return the value of the entry
   * \p (i,j).  This may be an
   * expensive operation, and you
   * should always be careful where
   * you call this function.
   */
  virtual T operator () (const libMesh::numeric_index_type i,
                         const libMesh::numeric_index_type j) const libmesh_override;

  /**
   * Return the l1-norm of the matrix, that is
   * \f$|M|_1=max_{all columns j}\sum_{all
   * rows i} |M_ij|\f$,
   * (max. sum of columns).
   * This is the
   * natural matrix norm that is compatible
   * to the l1-norm for vectors, i.e.
   * \f$|Mv|_1\leq |M|_1 |v|_1\f$.
   */
  virtual libMesh::Real l1_norm () const libmesh_override;

  /**
   * Return the linfty-norm of the
   * matrix, that is
   * \f$|M|_\infty=max_{all rows i}\sum_{all
   * columns j} |M_ij|\f$,
   * (max. sum of rows).
   * This is the
   * natural matrix norm that is compatible
   * to the linfty-norm of vectors, i.e.
   * \f$|Mv|_\infty \leq |M|_\infty |v|_\infty\f$.
   */
  virtual libMesh::Real linfty_norm () const libmesh_override;

  /**
   * see if GSL matrix has been closed
   * and fully assembled yet
   */
  virtual bool closed() const libmesh_override { return _closed; }

  /**
   * Print the contents of the matrix, by default to libMesh::out.
   * Currently identical to \p print().
   */
  virtual void print_personal(std::ostream & os=libMesh::out) const libmesh_override { this->print(os); }

  /**
   * Copies the diagonal part of the matrix into \p dest.
   */
  virtual void get_diagonal (libMesh::NumericVector<T> & dest) const libmesh_override;

  /**
   * Copies the transpose of the matrix into \p dest, which may be
   * *this.
   */
  virtual void get_transpose (libMesh::SparseMatrix<T> & dest) const libmesh_override;

private:

  std::unique_ptr<QUESO::BaseEnvironment> queso_env;
  QUESO::MpiComm queso_mpi_comm;
  std::unique_ptr<QUESO::Map> queso_map;

  /**
   * Actual QUESO::GslBlockMatrix we are wrapping.
   */
  DataType _mat;

  /**
   * Flag indicating if the matrix has been closed yet.
   */
  bool _closed;

  /**
   * Make other Eigen datatypes friends
   */
  friend class GslNumericVector<T>;
  // friend class EigenSparseLinearSolver<T>;
};

template <typename T>
GslNumericVector<T> operator*(const GslNumericBlockMatrix<T> & mat,
                              const GslNumericVector<T> & vec);

template <typename T>
GslNumericBlockMatrix<T> operator*(double a,
                             const GslNumericBlockMatrix<T> & mat);

template <typename T>
GslNumericBlockMatrix<T> matrixProduct(const GslNumericVector<T> & v1,
                                 const GslNumericVector<T> & v2);

template <typename T>
GslNumericBlockMatrix<T> operator*(const GslNumericBlockMatrix<T> & m1,
                             const GslNumericBlockMatrix<T> & m2);

template <typename T>
GslNumericBlockMatrix<T> operator+(const GslNumericBlockMatrix<T> & m1,
                             const GslNumericBlockMatrix<T> & m2);

} // namespace QUESO

#endif // #ifdef QUESO_GSL_NUMERIC_BLOCK_MATRIX_H