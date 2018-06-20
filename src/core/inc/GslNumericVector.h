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

#ifndef QUESO_GSL_NUMERIC_VECTOR
#define QUESO_GSL_NUMERIC_VECTOR

#include "libmesh/libmesh_common.h"
#include "libmesh/numeric_vector.h"

#include <queso/GslVector.h>

namespace libMesh
{

// Forward declarations
// template <typename T> class EigenSparseLinearSolver;
template <typename T> class SparseMatrix;

}

namespace QUESO
{

// Forward declarations
template <typename T> class GslSparseMatrix;

/**
 * Shim class for QUESO::GslVector
 */
template <typename T>
class GslNumericVector libmesh_final : public libMesh::NumericVector<T>
{
public:

  /**
   *  Dummy-Constructor. Dimension=0
   */
  explicit
  GslNumericVector (const libMesh::Parallel::Communicator & comm_in,
                     const libMesh::ParallelType = libMesh::AUTOMATIC);

  /**
   * Constructor. Set dimension to \p n and initialize all elements with zero.
   */
  explicit
  GslNumericVector (const libMesh::Parallel::Communicator & comm_in,
                     const libMesh::numeric_index_type n,
                     const libMesh::ParallelType = libMesh::AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global dimension
   * to \p n, and initialize all elements with zero.
   */
  GslNumericVector (const libMesh::Parallel::Communicator & comm_in,
                     const libMesh::numeric_index_type n,
                     const libMesh::numeric_index_type n_local,
                     const libMesh::ParallelType = libMesh::AUTOMATIC);

  /**
   * Constructor. Set local dimension to \p n_local, the global
   * dimension to \p n, but additionally reserve memory for the
   * indices specified by the \p ghost argument.
   */
  GslNumericVector (const libMesh::Parallel::Communicator & comm_in,
                     const libMesh::numeric_index_type N,
                     const libMesh::numeric_index_type n_local,
                     const std::vector<libMesh::numeric_index_type> & ghost,
                     const libMesh::ParallelType = libMesh::AUTOMATIC);

  /**
   * QUESO's ctors and methods
   */
  GslNumericVector(const BaseEnvironment & env, const Map & map);
  GslNumericVector(const GslNumericVector<T> & other);
  void cwSet(double value);
  void cwSet(unsigned int initialPos, const GslNumericVector<T> & vec);
  unsigned int sizeGlobal() const;
  unsigned int sizeLocal() const;
  const double & operator[](unsigned int i) const;
  double & operator[](unsigned int i);
  bool atLeastOneComponentSmallerThan(const GslNumericVector<T> & rhs) const;
  bool atLeastOneComponentBiggerThan(const GslNumericVector<T> & rhs) const;
  void cwExtract(unsigned int initialPos, GslNumericVector<T> & vec) const;
  double getMinValue() const;
  void cwSqrt();

  /**
   * Destructor, deallocates memory. Made virtual to allow
   * for derived classes to behave properly.
   */
  ~GslNumericVector ();

  /**
   * Convenient typedefs
   */
  typedef std::unique_ptr<GslVector> DataType;

  /**
   * Call the assemble functions
   */
  virtual void close () libmesh_override;

  /**
   * @returns the \p GslNumericVector to a pristine state.
   */
  virtual void clear () libmesh_override;

  /**
   * Set all entries to zero. Equivalent to \p v = 0, but more obvious and
   * faster.
   */
  virtual void zero () libmesh_override;

  /**
   * Creates a vector which has the same type, size and partitioning
   * as this vector, but whose data is all zero.  Returns it in an \p
   * UniquePtr.
   */
  virtual libMesh::UniquePtr<libMesh::NumericVector<T> > zero_clone () const libmesh_override;

  /**
   * Creates a copy of this vector and returns it in an \p UniquePtr.
   */
  virtual libMesh::UniquePtr<libMesh::NumericVector<T> > clone () const libmesh_override;

  /**
   * Change the dimension of the vector to \p N. The reserved memory for
   * this vector remains unchanged if possible, to make things faster, but
   * this may waste some memory, so take this in the back of your head.
   * However, if \p N==0 all memory is freed, i.e. if you want to resize
   * the vector and release the memory not needed, you have to first call
   * \p init(0) and then \p init(N). This cited behaviour is analogous
   * to that of the STL containers.
   *
   * On \p fast==false, the vector is filled by
   * zeros.
   */
  virtual void init (const libMesh::numeric_index_type N,
                     const libMesh::numeric_index_type n_local,
                     const bool         fast=false,
                     const libMesh::ParallelType ptype=libMesh::AUTOMATIC) libmesh_override;

  /**
   * call init with n_local = N,
   */
  virtual void init (const libMesh::numeric_index_type N,
                     const bool         fast=false,
                     const libMesh::ParallelType ptype=libMesh::AUTOMATIC) libmesh_override;

  /**
   * Create a vector that holds tha local indices plus those specified
   * in the \p ghost argument.
   */
  virtual void init (const libMesh::numeric_index_type /*N*/,
                     const libMesh::numeric_index_type /*n_local*/,
                     const std::vector<libMesh::numeric_index_type> & /*ghost*/,
                     const bool /*fast*/ = false,
                     const libMesh::ParallelType = libMesh::AUTOMATIC) libmesh_override;

  /**
   * Creates a vector that has the same dimension and storage type as
   * \p other, including ghost dofs.
   */
  virtual void init (const libMesh::NumericVector<T> & other,
                     const bool fast = false) libmesh_override;

  /**
   * \f$U(0-N) = s\f$: fill all components.
   */
  virtual libMesh::NumericVector<T> & operator= (const T s) libmesh_override;

  /**
   *  \f$U = V\f$: copy all components.
   */
  virtual libMesh::NumericVector<T> & operator= (const libMesh::NumericVector<T> & v) libmesh_override;

  /**
   *  \f$U = V\f$: copy all components.
   */
  GslNumericVector<T> & operator= (const GslNumericVector<T> & v);

  /**
   *  \f$U = V\f$: copy all components.
   */
  virtual libMesh::NumericVector<T> & operator= (const std::vector<T> & v) libmesh_override;

  /**
   * @returns the minimum element in the vector.
   * In case of complex numbers, this returns the minimum
   * Real part.
   */
  virtual libMesh::Real min () const libmesh_override;

  /**
   * @returns the maximum element in the vector.
   * In case of complex numbers, this returns the maximum
   * Real part.
   */
  virtual libMesh::Real max () const libmesh_override;

  /**
   * @returns the sum of values in a vector
   */
  virtual T sum () const libmesh_override;

  /**
   * @returns the \f$l_1\f$-norm of the vector, i.e.
   * the sum of the absolute values.
   */
  virtual libMesh::Real l1_norm () const libmesh_override;

  /**
   * @returns the \f$l_2\f$-norm of the vector, i.e.
   * the square root of the sum of the
   * squares of the elements.
   */
  virtual libMesh::Real l2_norm () const libmesh_override;

  /**
   * @returns the maximum absolute value of the
   * elements of this vector, which is the
   * \f$l_\infty\f$-norm of a vector.
   */
  virtual libMesh::Real linfty_norm () const libmesh_override;

  /**
   * @returns dimension of the vector. This
   * function was formerly called \p n(), but
   * was renamed to get the \p GslNumericVector class
   * closer to the C++ standard library's
   * \p std::vector container.
   */
  virtual libMesh::numeric_index_type size () const libmesh_override;

  /**
   * @returns the local size of the vector
   * (index_stop-index_start)
   */
  virtual libMesh::numeric_index_type local_size() const libmesh_override;

  /**
   * @returns the index of the first vector element
   * actually stored on this processor
   */
  virtual libMesh::numeric_index_type first_local_index() const libmesh_override;

  /**
   * @returns the index of the last vector element
   * actually stored on this processor
   */
  virtual libMesh::numeric_index_type last_local_index() const libmesh_override;

  /**
   * Access components, returns \p U(i).
   */
  virtual T operator() (const libMesh::numeric_index_type i) const libmesh_override;

  /**
   * Addition operator.
   * Fast equivalent to \p U.add(1, V).
   */
  virtual libMesh::NumericVector<T> & operator += (const libMesh::NumericVector<T> & v) libmesh_override;

  /**
   * Subtraction operator.
   * Fast equivalent to \p U.add(-1, V).
   */
  virtual libMesh::NumericVector<T> & operator -= (const libMesh::NumericVector<T> & v) libmesh_override;

  /**
   * Pointwise Division operator. ie divide every entry in this vector by the entry in v
   */
  virtual libMesh::NumericVector<T> & operator /= (const libMesh::NumericVector<T> & v_in) libmesh_override;

  /**
   * Replace each entry v_i of this vector by its reciprocal, 1/v_i.
   */
  virtual void reciprocal() libmesh_override;

  /**
   * Replace each entry v_i = real(v_i) + imag(v_i)
   * of this vector by its complex conjugate, real(v_i) - imag(v_i)
   */
  virtual void conjugate() libmesh_override;

  /**
   * v(i) = value
   */
  virtual void set (const libMesh::numeric_index_type i, const T value) libmesh_override;

  /**
   * v(i) += value
   */
  virtual void add (const libMesh::numeric_index_type i, const T value) libmesh_override;

  /**
   * \f$U(0-LIBMESH_DIM)+=s\f$.
   * Addition of \p s to all components. Note
   * that \p s is a scalar and not a vector.
   */
  virtual void add (const T s) libmesh_override;

  /**
   * \f$ U+=V \f$.
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  virtual void add (const libMesh::NumericVector<T> & V) libmesh_override;

  /**
   * \f$ U+=a*V \f$.
   * Simple vector addition, equal to the
   * \p operator +=.
   */
  virtual void add (const T a, const libMesh::NumericVector<T> & v) libmesh_override;

  /**
   * We override one NumericVector<T>::add_vector() method but don't
   * want to hide the other defaults.
   */
  using libMesh::NumericVector<T>::add_vector;

  /**
   * \f$U+=A*V\f$, add the product of a \p SparseMatrix \p A
   * and a \p NumericVector \p V to \p this, where \p this=U.
   */
  virtual void add_vector (const libMesh::NumericVector<T> &,
                           const libMesh::SparseMatrix<T> &) libmesh_override;

  /**
   * \f$U+=A^T*V\f$, add the product of the transpose of a \p SparseMatrix \p A_trans
   * and a \p NumericVector \p V to \p this, where \p this=U.
   */
  virtual void add_vector_transpose (const libMesh::NumericVector<T> &,
                                     const libMesh::SparseMatrix<T> &) libmesh_override;

  /**
   * Scale each element of the
   * vector by the given factor.
   */
  virtual void scale (const T factor) libmesh_override;

  /**
   * v = abs(v)... that is, each entry in v is replaced
   * by its absolute value.
   */
  virtual void abs() libmesh_override;

  /**
   * Computes the dot product, p = U.V
   */
  virtual T dot(const libMesh::NumericVector<T> & V) const libmesh_override;

  /**
   * Creates a copy of the global vector in the
   * local vector \p v_local.
   */
  virtual void localize (std::vector<T> & v_local) const libmesh_override;

  /**
   * Same, but fills a \p NumericVector<T> instead of
   * a \p std::vector.
   */
  virtual void localize (libMesh::NumericVector<T> & v_local) const libmesh_override;

  /**
   * Creates a local vector \p v_local containing
   * only information relevant to this processor, as
   * defined by the \p send_list.
   */
  virtual void localize (libMesh::NumericVector<T> & v_local,
                         const std::vector<libMesh::numeric_index_type> & send_list) const libmesh_override;

  /**
   * Fill in the local std::vector "v_local" with the global indices
   * given in "indices".  See numeric_vector.h for more details.
   */
  virtual void localize (std::vector<T> & v_local,
                         const std::vector<libMesh::numeric_index_type> & indices) const libmesh_override;

  /**
   * Updates a local vector with selected values from neighboring
   * processors, as defined by \p send_list.
   */
  virtual void localize (const libMesh::numeric_index_type first_local_idx,
                         const libMesh::numeric_index_type last_local_idx,
                         const std::vector<libMesh::numeric_index_type> & send_list) libmesh_override;

  /**
   * Creates a local copy of the global vector in
   * \p v_local only on processor \p proc_id.  By
   * default the data is sent to processor 0.  This method
   * is useful for outputting data from one processor.
   */
  virtual void localize_to_one (std::vector<T> & v_local,
                                const libMesh::processor_id_type proc_id=0) const libmesh_override;

  /**
   * Computes the pointwise (i.e. component-wise) product of \p vec1
   * and \p vec2 and stores the result in \p *this.
   */
  virtual void pointwise_mult (const libMesh::NumericVector<T> & vec1,
                               const libMesh::NumericVector<T> & vec2) libmesh_override;

  /**
   * Swaps the contents.
   */
  virtual void swap (libMesh::NumericVector<T> & v) libmesh_override;

  /**
   * References to the underlying QUESO-wrapped GSL data types. Note this is generally
   * not required in user-level code.
   */
  QUESO::GslVector &       vec ()        { return *_vec; }
  const QUESO::GslVector & vec () const  { return *_vec; }

private:

  static std::map<const QUESO::MpiComm *, libMesh::Parallel::Communicator> comm_map;

  std::unique_ptr<QUESO::BaseEnvironment> queso_env;
  QUESO::MpiComm queso_mpi_comm;
  std::unique_ptr<QUESO::Map> queso_map;

  /**
   * Actual QUESO::GslVector we are wrapping.
   */
  DataType _vec;

  /**
   * Make other QUESO-wrapped GSL datatypes friends
   */
  friend class GslSparseMatrix<T>;
  // friend class libMesh::EigenSparseLinearSolver<T>;
};



// ---------------------------------------------------------
// GslNumericVector inline methods
template <typename T>
inline
GslNumericVector<T>::GslNumericVector (const libMesh::Parallel::Communicator & comm_in,
                                         const libMesh::ParallelType ptype)
  : libMesh::NumericVector<T>(comm_in, ptype),
    queso_env(new EmptyEnvironment()),
    queso_mpi_comm(*queso_env, comm_in.get())
{
  this->_type = ptype;
}



template <typename T>
inline
GslNumericVector<T>::GslNumericVector (const libMesh::Parallel::Communicator & comm_in,
                                         const libMesh::numeric_index_type n,
                                         const libMesh::ParallelType ptype)
  : libMesh::NumericVector<T>(comm_in, ptype),
    queso_env(new EmptyEnvironment()),
    queso_mpi_comm(*queso_env, comm_in.get())
{
  this->init(n, n, false, ptype);
}



template <typename T>
inline
GslNumericVector<T>::GslNumericVector (const libMesh::Parallel::Communicator & comm_in,
                                         const libMesh::numeric_index_type n,
                                         const libMesh::numeric_index_type n_local,
                                         const libMesh::ParallelType ptype)
  : libMesh::NumericVector<T>(comm_in, ptype),
    queso_env(new EmptyEnvironment()),
    queso_mpi_comm(*queso_env, comm_in.get())
{
  this->init(n, n_local, false, ptype);
}



template <typename T>
inline
GslNumericVector<T>::GslNumericVector (const libMesh::Parallel::Communicator & comm_in,
                                         const libMesh::numeric_index_type N,
                                         const libMesh::numeric_index_type n_local,
                                         const std::vector<libMesh::numeric_index_type> & ghost,
                                         const libMesh::ParallelType ptype)
  : libMesh::NumericVector<T>(comm_in, ptype),
    queso_env(new EmptyEnvironment()),
    queso_mpi_comm(*queso_env, comm_in.get())
{
  this->init(N, n_local, ghost, false, ptype);
}



template <typename T>
inline
GslNumericVector<T>::~GslNumericVector ()
{
  this->clear ();
}



template <typename T>
inline
void GslNumericVector<T>::init (const libMesh::numeric_index_type n,
                                 const libMesh::numeric_index_type n_local,
                                 const bool fast,
                                 const libMesh::ParallelType)
{
  // GSL vectors only for serial cases,
  // but can provide a "parallel" vector on one processor.
  if (n != n_local)
    libmesh_error_msg("Error: GslNumericVectors can only be used in serial!");

  this->_type = libMesh::SERIAL;

  // Clear initialized vectors
  if (this->initialized())
    this->clear();

  this->queso_map.reset(new QUESO::Map(n, 0, this->queso_mpi_comm));
  this->_vec.reset(new QUESO::GslVector(*(this->queso_env), *(this->queso_map)));

  this->_is_initialized = true;
#ifndef NDEBUG
  this->_is_closed = true;
#endif

  // Optionally zero out all components
  if (fast == false)
    this->zero ();

  return;
}



template <typename T>
inline
void GslNumericVector<T>::init (const libMesh::numeric_index_type n,
                                 const bool fast,
                                 const libMesh::ParallelType ptype)
{
  this->init(n,n,fast,ptype);
}


template <typename T>
inline
void GslNumericVector<T>::init (const libMesh::numeric_index_type n,
                                 const libMesh::numeric_index_type n_local,
                                 const std::vector<libMesh::numeric_index_type> & libmesh_dbg_var(ghost),
                                 const bool fast,
                                 const libMesh::ParallelType ptype)
{
  libmesh_assert(ghost.empty());
  this->init(n,n_local,fast,ptype);
}



/* Default implementation for solver packages for which ghosted
   vectors are not yet implemented.  */
template <class T>
void GslNumericVector<T>::init (const libMesh::NumericVector<T> & other,
                                 const bool fast)
{
  this->init(other.size(),other.local_size(),fast,other.type());
}



template <typename T>
inline
void GslNumericVector<T>::close ()
{
  libmesh_assert (this->initialized());

#ifndef NDEBUG
  this->_is_closed = true;
#endif
}



template <typename T>
inline
void GslNumericVector<T>::clear ()
{
  this->queso_map.reset(new QUESO::Map(0, 0, this->queso_mpi_comm));
  this->_vec.reset(new QUESO::GslVector(*(this->queso_env), *(this->queso_map)));

  this->_is_initialized = false;
#ifndef NDEBUG
  this->_is_closed = false;
#endif
}



template <typename T> inline
void GslNumericVector<T>::zero ()
{
  libmesh_assert (this->initialized());
  libmesh_assert (this->closed());

  _vec->cwSet(0.0);
}



template <typename T>
inline
libMesh::UniquePtr<libMesh::NumericVector<T> > GslNumericVector<T>::zero_clone () const
{
  libMesh::NumericVector<T> * cloned_vector = new GslNumericVector<T>(this->comm());
  cloned_vector->init(*this);
  return libMesh::UniquePtr<libMesh::NumericVector<T> >(cloned_vector);
}



template <typename T>
inline
libMesh::UniquePtr<libMesh::NumericVector<T> > GslNumericVector<T>::clone () const
{
  libMesh::NumericVector<T> * cloned_vector = new GslNumericVector<T>(this->comm());
  cloned_vector->init(*this, true);
  *cloned_vector = *this;
  return libMesh::UniquePtr<libMesh::NumericVector<T> >(cloned_vector);
}



template <typename T>
inline
libMesh::numeric_index_type GslNumericVector<T>::size () const
{
  libmesh_assert (this->initialized());

  return static_cast<libMesh::numeric_index_type>(_vec->sizeLocal());
}



template <typename T>
inline
libMesh::numeric_index_type GslNumericVector<T>::local_size () const
{
  libmesh_assert (this->initialized());

  return this->size();
}



template <typename T>
inline
libMesh::numeric_index_type GslNumericVector<T>::first_local_index () const
{
  libmesh_assert (this->initialized());

  return 0;
}



template <typename T>
inline
libMesh::numeric_index_type GslNumericVector<T>::last_local_index () const
{
  libmesh_assert (this->initialized());

  return this->size();
}



template <typename T>
inline
void GslNumericVector<T>::set (const libMesh::numeric_index_type i, const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->size());

  (*_vec)[static_cast<int>(i)] = value;

#ifndef NDEBUG
  this->_is_closed = false;
#endif
}



template <typename T>
inline
void GslNumericVector<T>::add (const libMesh::numeric_index_type i, const T value)
{
  libmesh_assert (this->initialized());
  libmesh_assert_less (i, this->size());

  (*_vec)[static_cast<int>(i)] += value;

#ifndef NDEBUG
  this->_is_closed = false;
#endif
}



template <typename T>
inline
T GslNumericVector<T>::operator() (const libMesh::numeric_index_type i) const
{
  libmesh_assert (this->initialized());
  libmesh_assert ( ((i >= this->first_local_index()) &&
                    (i <  this->last_local_index())) );

  return (*_vec)[static_cast<int>(i)];
}



template <typename T>
inline
void GslNumericVector<T>::swap (libMesh::NumericVector<T> & other)
{
  GslNumericVector<T> & v = libMesh::cast_ref<GslNumericVector<T> &>(other);

  // Copy other to *tmp
  GslNumericVector<T> * tmp = libMesh::cast_ptr<GslNumericVector<T> *>(other.clone().get());

  // Store *this in other
  *v._vec = *_vec;

  // Copy *tmp to *this
  *_vec = *(tmp->_vec);

  std::swap (this->_is_closed,      v._is_closed);
  std::swap (this->_is_initialized, v._is_initialized);
  std::swap (this->_type,           v._type);
}


} // namespace libMesh


#endif  // QUESO_GSL_NUMERIC_VECTOR
