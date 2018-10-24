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
#include <algorithm> // for std::min
#include <limits>

// Local Includes
#include "libmesh/dense_subvector.h"
#include "libmesh/dense_vector.h"
#include "libmesh/eigen_sparse_vector.h"
#include "libmesh/eigen_sparse_matrix.h"

// Just so my fingers don't hate me after editing 300+ source files
using namespace QUESO;

#ifdef LIBMESH_HAVE_EIGEN

namespace libMesh
{

template <typename T>
EigenSparseVector<T>::EigenSparseVector(const BaseEnvironment & env,
                                        const Map & map)
  : libMesh::NumericVector<T>(
      comm_map.emplace(std::make_pair(
          &(map.Comm()),
          libMesh::Parallel::Communicator(map.Comm().Comm()))).first->second,
      libMesh::AUTOMATIC),
    queso_env(new EmptyEnvironment()),
    queso_mpi_comm(map.Comm())
{
  this->init(map.NumGlobalElements(),
             map.NumGlobalElements(),
             false,
             libMesh::AUTOMATIC);
}

template <typename T>
EigenSparseVector<T>::EigenSparseVector(const EigenSparseVector<T> & other)
  : libMesh::NumericVector<T>(
      comm_map.emplace(std::make_pair(
          &(other.queso_map->Comm()),
          libMesh::Parallel::Communicator(other.queso_map->Comm().Comm()))).first->second,
      libMesh::AUTOMATIC),
    queso_env(new EmptyEnvironment()),  // This is empty but we don't care
                                        // because only the internal queso data
                                        // type cares about the environment
    queso_mpi_comm(other.queso_map->Comm())
{
  this->init(other.queso_map->NumGlobalElements(),
             other.queso_map->NumGlobalElements(),
             false,
             libMesh::AUTOMATIC);
}

template <typename T>
void
EigenSparseVector<T>::cwSet(double value)
{
  libmesh_assert (this->initialized());
  libmesh_assert (this->closed());

  _vec.fill(value);
}

template <typename T>
unsigned int
EigenSparseVector<T>::sizeGlobal() const
{
  return this->size();
}

template <typename T>
unsigned int
EigenSparseVector<T>::sizeLocal() const
{
  return this->size();
}

template <typename T>
const double &
EigenSparseVector<T>::operator[](unsigned int i) const
{
  return this->_vec[static_cast<eigen_idx_type>(i)];
}

template <typename T>
double &
EigenSparseVector<T>::operator[](unsigned int i)
{
  return this->_vec[static_cast<eigen_idx_type>(i)];
}

template <typename T>
bool
EigenSparseVector<T>::atLeastOneComponentSmallerThan(const EigenSparseVector<T> & rhs) const
{
  queso_require_equal_to_msg(this->sizeLocal(), rhs.sizeLocal(), "vectors have different sizes");

  bool result = false;
  unsigned int i = 0;
  unsigned int size = this->sizeLocal();
  while ((i < size) && (result == false)) {
    result = ( (*this)[i] < rhs[i] );
    i++;
  };

  return result;
}

template <typename T>
bool
EigenSparseVector<T>::atLeastOneComponentBiggerThan(const EigenSparseVector<T> & rhs) const
{
  queso_require_equal_to_msg(this->sizeLocal(), rhs.sizeLocal(), "vectors have different sizes");

  bool result = false;
  unsigned int i = 0;
  unsigned int size = this->sizeLocal();
  while ((i < size) && (result == false)) {
    result = ( (*this)[i] > rhs[i] );
    i++;
  };

  return result;
}

template <typename T>
void
EigenSparseVector<T>::cwExtract(unsigned int initialPos, EigenSparseVector<T> & vec) const
{
  queso_require_less_msg(initialPos, this->sizeLocal(), "invalid initialPos");
  queso_require_less_equal_msg((initialPos +vec.sizeLocal()), this->sizeLocal(), "invalid vec.sizeLocal()");

  for (unsigned int i = 0; i < vec.sizeLocal(); ++i) {
    vec[i] = (*this)[initialPos+i];
  }
}

template <typename T>
void
EigenSparseVector<T>::cwSet(unsigned int initialPos, const EigenSparseVector<T> & vec)
{
  queso_require_less_msg(initialPos, this->sizeLocal(), "invalid initialPos");
  queso_require_less_equal_msg((initialPos +vec.sizeLocal()), this->sizeLocal(), "invalid vec.sizeLocal()");

  for (unsigned int i = 0; i < vec.sizeLocal(); ++i) {
    (*this)[initialPos+i] = vec[i];
  }
}

template <typename T>
void
EigenSparseVector<T>::cwSetGaussian(const EigenSparseVector<T> & meanVec, const EigenSparseVector<T> & stdDevVec)
{
  // should be a non-member non-friend function
  queso_not_implemented();
}

template <typename T>
void
EigenSparseVector<T>::cwSetUniform(const EigenSparseVector<T> & aVec, const EigenSparseVector<T> & bVec)
{
  // should be a non-member non-friend function
  queso_not_implemented();
}

template <typename T>
double
EigenSparseVector<T>::getMinValue( ) const
{
  return _vec.minCoeff();
}

template <typename T>
void
EigenSparseVector<T>::setPrintHorizontally(bool value) const
{
  m_printHorizontally = value;
}

template <typename T>
bool
EigenSparseVector<T>::getPrintHorizontally() const
{
  return m_printHorizontally;
}

template <typename T>
void
EigenSparseVector<T>::setPrintScientific(bool value) const
{
  m_printScientific = value;
}

template <typename T>
bool
EigenSparseVector<T>::getPrintScientific() const
{
  return m_printScientific;
}

template <typename T>
int
EigenSparseVector<T>::getMaxValueIndex() const
{
  T maxval = -INFINITY;  // Assumes T is floating point type?
  unsigned int maxval_index = 0;

  for (unsigned int i = 0; i < this->size(); i++) {
    if (_vec[i] > maxval) {
      maxval = _vec[i];
      maxval_index = i;
    }
  }

  return maxval_index;
}

template <typename T>
double
EigenSparseVector<T>::norm2() const
{
  return l2_norm();
}

template <typename T>
unsigned int
EigenSparseVector<T>::numOfProcsForStorage() const
{
  return queso_map->Comm().NumProc();
}

template <typename T>
void
EigenSparseVector<T>::cwSetBeta(const EigenSparseVector<T> & a, const EigenSparseVector<T> & b)
{
  // should be a non-member non-friend function
  queso_not_implemented();
}

template <typename T>
void
EigenSparseVector<T>::cwSetConcatenated(const std::vector<const EigenSparseVector<T>* >& vecs)
{
  unsigned int cummulativeSize = 0;
  for (unsigned int i = 0; i < vecs.size(); ++i) {
    EigenSparseVector<T> tmpVec(*(vecs[i]));
    for (unsigned int j = 0; j < vecs[i]->sizeLocal(); ++j) {
      (*this)[cummulativeSize+j] = tmpVec[j];
    }
    cummulativeSize += vecs[i]->sizeLocal();
  }

  queso_require_equal_to_msg(this->sizeLocal(), cummulativeSize, "incompatible vector sizes");
}

template <typename T>
void
EigenSparseVector<T>::cwSetGamma(const EigenSparseVector<T> & a, const EigenSparseVector<T> & b)
{
  queso_not_implemented();
}

template <typename T>
void
EigenSparseVector<T>::cwSetGaussian(double mean, double stdDev)
{
  queso_not_implemented();
}

template <typename T>
void
EigenSparseVector<T>::cwSetInverseGamma(const EigenSparseVector<T> & a, const EigenSparseVector<T> & b)
{
  queso_not_implemented();
}

template <typename T>
double
EigenSparseVector<T>::sumOfComponents() const
{
  return _vec.sum();
}

template <typename T>
bool
EigenSparseVector<T>::atLeastOneComponentSmallerOrEqualThan(const EigenSparseVector<T> & rhs) const
{
  queso_require_equal_to_msg(this->sizeLocal(), rhs.sizeLocal(), "vectors have different sizes");

  bool result = false;
  unsigned int i = 0;
  unsigned int size = this->sizeLocal();
  while ((i < size) && (result == false)) {
    result = ( (*this)[i] <= rhs[i] );
    i++;
  };

  return result;
}

template <typename T>
T EigenSparseVector<T>::sum () const
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());

  return _vec.sum();
}



template <typename T>
Real EigenSparseVector<T>::l1_norm () const
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());

  return _vec.lpNorm<1>();
}



template <typename T>
Real EigenSparseVector<T>::l2_norm () const
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());

  return _vec.lpNorm<2>();
}



template <typename T>
Real EigenSparseVector<T>::linfty_norm () const
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());

  return _vec.lpNorm<Eigen::Infinity>();
}



template <typename T>
NumericVector<T> & EigenSparseVector<T>::operator += (const NumericVector<T> & v_in)
{
  libmesh_assert (this->closed());

  const EigenSparseVector<T> & v = cast_ref<const EigenSparseVector<T> &>(v_in);

  _vec += v._vec;

  return *this;
}




template <typename T>
NumericVector<T> & EigenSparseVector<T>::operator -= (const NumericVector<T> & v_in)
{
  libmesh_assert (this->closed());

  const EigenSparseVector<T> & v = cast_ref<const EigenSparseVector<T> &>(v_in);

  _vec -= v._vec;

  return *this;
}



template <typename T>
NumericVector<T> & EigenSparseVector<T>::operator /= (const NumericVector<T> & v_in)
{
  libmesh_assert (this->closed());
  libmesh_assert_equal_to(size(), v_in.size());

  const EigenSparseVector<T> & v = cast_ref<const EigenSparseVector<T> &>(v_in);

  _vec = _vec.cwiseQuotient(v._vec);

  return *this;
}




template <typename T>
void EigenSparseVector<T>::reciprocal()
{
#ifndef NDEBUG
  const numeric_index_type n = this->size();

  for (numeric_index_type i=0; i<n; i++)
    // Don't divide by zero!
    libmesh_assert_not_equal_to ((*this)(i), T(0));
#endif

  _vec = _vec.cwiseInverse();
}



template <typename T>
void EigenSparseVector<T>::conjugate()
{
  _vec = _vec.conjugate();
}



template <typename T>
void EigenSparseVector<T>::add (const T v)
{
  _vec += EigenSV::Constant(this->size(), v);

#ifndef NDEBUG
  this->_is_closed = false;
#endif
}




template <typename T>
void EigenSparseVector<T>::add (const NumericVector<T> & v_in)
{
  libmesh_assert (this->initialized());

  const EigenSparseVector<T> & v = cast_ref<const EigenSparseVector<T> &>(v_in);

  _vec += v._vec;
}



template <typename T>
void EigenSparseVector<T>::add (const T a, const NumericVector<T> & v_in)
{
  libmesh_assert (this->initialized());

  const EigenSparseVector<T> & v = cast_ref<const EigenSparseVector<T> &>(v_in);

  _vec += v._vec*a;
}



template <typename T>
void EigenSparseVector<T>::add_vector (const NumericVector<T> & vec_in,
                                       const SparseMatrix<T>  & mat_in)
{
  // Make sure the data passed in are really in Eigen types
  const EigenSparseVector<T> * e_vec = cast_ptr<const EigenSparseVector<T> *>(&vec_in);
  const EigenSparseMatrix<T> * mat = cast_ptr<const EigenSparseMatrix<T> *>(&mat_in);

  libmesh_assert(e_vec);
  libmesh_assert(mat);

  _vec += mat->_mat*e_vec->_vec;
}



template <typename T>
void EigenSparseVector<T>::add_vector_transpose (const NumericVector<T> & vec_in,
                                                 const SparseMatrix<T>  & mat_in)
{
  // Make sure the data passed in are really in Eigen types
  const EigenSparseVector<T> * e_vec = cast_ptr<const EigenSparseVector<T> *>(&vec_in);
  const EigenSparseMatrix<T> * mat = cast_ptr<const EigenSparseMatrix<T> *>(&mat_in);

  libmesh_assert(e_vec);
  libmesh_assert(mat);

  _vec += mat->_mat.transpose()*e_vec->_vec;
}



template <typename T>
void EigenSparseVector<T>::scale (const T factor)
{
  libmesh_assert (this->initialized());

  _vec *= factor;
}



template <typename T>
void EigenSparseVector<T>::abs()
{
  libmesh_assert (this->initialized());

  const numeric_index_type n = this->size();

  for (numeric_index_type i=0; i!=n; ++i)
    this->set(i,std::abs((*this)(i)));
}



template <typename T>
T EigenSparseVector<T>::dot (const NumericVector<T> & V) const
{
  libmesh_assert (this->initialized());

  // Make sure the NumericVector passed in is really a EigenSparseVector
  const EigenSparseVector<T> * v = cast_ptr<const EigenSparseVector<T> *>(&V);
  libmesh_assert(v);

  return _vec.dot(v->_vec);
}



template <typename T>
NumericVector<T> &
EigenSparseVector<T>::operator = (const T s)
{
  libmesh_assert (this->initialized());
  libmesh_assert (this->closed());

  _vec.fill(s);

  return *this;
}



template <typename T>
NumericVector<T> &
EigenSparseVector<T>::operator = (const NumericVector<T> & v_in)
{
  // Make sure the NumericVector passed in is really a EigenSparseVector
  const EigenSparseVector<T> * v =
    cast_ptr<const EigenSparseVector<T> *>(&v_in);

  libmesh_assert(v);

  *this = *v;

  return *this;
}



template <typename T>
EigenSparseVector<T> &
EigenSparseVector<T>::operator = (const EigenSparseVector<T> & v)
{
  libmesh_assert (this->initialized());
  libmesh_assert (v.closed());
  libmesh_assert_equal_to (this->size(), v.size());

  _vec = v._vec;

#ifndef NDEBUG
  this->_is_closed = true;
#endif

  return *this;
}



template <typename T>
NumericVector<T> &
EigenSparseVector<T>::operator = (const std::vector<T> & v)
{
  /**
   * Case 1:  The vector is the same size of
   * The global vector.  Only add the local components.
   */
  if (this->size() == v.size())
    for (numeric_index_type i=0; i<v.size(); i++)
      this->set (i, v[i]);

  else
    libmesh_error_msg("this->size() = " << this->size() << " must be equal to v.size() = " << v.size());

  return *this;
}


template <typename T>
void EigenSparseVector<T>::localize (NumericVector<T> & v_local_in) const
{
  // Make sure the NumericVector passed in is really a EigenSparseVector
  EigenSparseVector<T> * v_local =
    cast_ptr<EigenSparseVector<T> *>(&v_local_in);

  libmesh_assert(v_local);

  *v_local = *this;
}



template <typename T>
void EigenSparseVector<T>::localize (NumericVector<T> & v_local_in,
                                     const std::vector<numeric_index_type> & libmesh_dbg_var(send_list)) const
{
  // Make sure the NumericVector passed in is really a EigenSparseVector
  EigenSparseVector<T> * v_local =
    cast_ptr<EigenSparseVector<T> *>(&v_local_in);

  libmesh_assert(v_local);
  libmesh_assert_less_equal (send_list.size(), v_local->size());

  *v_local = *this;
}



template <typename T>
void EigenSparseVector<T>::localize (std::vector<T> & v_local,
                                     const std::vector<numeric_index_type> & indices) const
{
  // EigenSparseVectors are serial, so we can just copy values
  v_local.resize(indices.size());

  for (numeric_index_type i=0; i<v_local.size(); i++)
    v_local[i] = (*this)(indices[i]);
}



template <typename T>
void EigenSparseVector<T>::localize (const numeric_index_type libmesh_dbg_var(first_local_idx),
                                     const numeric_index_type libmesh_dbg_var(last_local_idx),
                                     const std::vector<numeric_index_type> & libmesh_dbg_var(send_list))
{
  libmesh_assert_equal_to (first_local_idx, 0);
  libmesh_assert_equal_to (last_local_idx+1, this->size());

  libmesh_assert_less_equal (send_list.size(), this->size());

#ifndef NDEBUG
  this->_is_closed = true;
#endif
}



template <typename T>
void EigenSparseVector<T>::localize (std::vector<T> & v_local) const

{
  v_local.resize(this->size());

  for (numeric_index_type i=0; i<v_local.size(); i++)
    v_local[i] = (*this)(i);
}



template <typename T>
void EigenSparseVector<T>::localize_to_one (std::vector<T> & v_local,
                                            const processor_id_type libmesh_dbg_var(pid)) const
{
  libmesh_assert_equal_to (pid, 0);

  this->localize (v_local);
}



template <typename T>
void EigenSparseVector<T>::pointwise_mult (const NumericVector<T> & /*vec1*/,
                                           const NumericVector<T> & /*vec2*/)
{
  libmesh_not_implemented();
}



template <typename T>
Real EigenSparseVector<T>::max() const
{
  libmesh_assert (this->initialized());
  if (!this->size())
    return -std::numeric_limits<Real>::max();

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  Real the_max = libmesh_real((*this)(0));

  const numeric_index_type n = this->size();

  for (numeric_index_type i=1; i<n; i++)
    the_max = std::max (the_max, libmesh_real((*this)(i)));

  return the_max;
#else
  return libmesh_real(_vec.maxCoeff());
#endif
}



template <typename T>
Real EigenSparseVector<T>::min () const
{
  libmesh_assert (this->initialized());
  if (!this->size())
    return std::numeric_limits<Real>::max();

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  Real the_min = libmesh_real((*this)(0));

  const numeric_index_type n = this->size();

  for (numeric_index_type i=1; i<n; i++)
    the_min = std::min (the_min, libmesh_real((*this)(i)));

  return the_min;
#else
  return libmesh_real(_vec.minCoeff());
#endif
}


//------------------------------------------------------------------
// Explicit instantiations
template class EigenSparseVector<Number>;

template <typename T>
std::map<const QUESO::MpiComm *, Parallel::Communicator>
EigenSparseVector<T>::comm_map;

} // namespace libMesh

// Not really sure what namespace I should be extending
namespace QUESO
{

template <typename T>
libMesh::EigenSparseVector<T>
operator+(const libMesh::EigenSparseVector<T> & x, const libMesh::EigenSparseVector<T> & y)
{
  libMesh::EigenSparseVector<T> answer(x);
  answer.vec() += y.vec();
  return answer;
}

template <typename T>
libMesh::EigenSparseVector<T> operator-(const libMesh::EigenSparseVector<T> & x,
                                        const libMesh::EigenSparseVector<T> & y)
{
  // How to write inefficient code lesson 1
  libMesh::EigenSparseVector<T> answer(x);
  answer -= y;
  return answer;
}

template <typename T>
libMesh::EigenSparseVector<T> operator*(double a,
                                        const libMesh::EigenSparseVector<T> & y)
{
  libMesh::EigenSparseVector<T> answer(y);
  answer.scale(a);
  return answer;
}

template <typename T>
libMesh::EigenSparseVector<T> operator*(const libMesh::EigenSparseVector<T> & x,
                                        const libMesh::EigenSparseVector<T> & y)
{
  libMesh::EigenSparseVector<T> answer(x);
  answer.vec() *= y.vec();
  return answer;
}

template <typename T>
libMesh::EigenSparseVector<T>
operator/(const libMesh::EigenSparseVector<T> & x, const libMesh::EigenSparseVector<T> & y)
{
  libMesh::EigenSparseVector<T> answer(x);
  answer.vec() = answer.vec().cwiseQuotient(y.vec());
  return answer;
}

template libMesh::EigenSparseVector<libMesh::Number> operator+(const libMesh::EigenSparseVector<libMesh::Number> & x, const libMesh::EigenSparseVector<libMesh::Number> & y);
template libMesh::EigenSparseVector<libMesh::Number> operator-(const libMesh::EigenSparseVector<libMesh::Number> & x, const libMesh::EigenSparseVector<libMesh::Number> & y);
template libMesh::EigenSparseVector<libMesh::Number> operator*(double a, const libMesh::EigenSparseVector<libMesh::Number> & y);
template libMesh::EigenSparseVector<libMesh::Number> operator*(const libMesh::EigenSparseVector<libMesh::Number> & x, const libMesh::EigenSparseVector<libMesh::Number> & y);
template libMesh::EigenSparseVector<libMesh::Number> operator/(const libMesh::EigenSparseVector<libMesh::Number> & x, const libMesh::EigenSparseVector<libMesh::Number> & y);

}


#endif // #ifdef LIBMESH_HAVE_EIGEN
