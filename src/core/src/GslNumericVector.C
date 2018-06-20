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

#include <queso/GslNumericVector.h>
#include <queso/GslSparseMatrix.h>

namespace QUESO
{

template <typename T>
T GslNumericVector<T>::sum () const
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());

  return _vec->sumOfComponents();
}



template <typename T>
libMesh::Real GslNumericVector<T>::l1_norm () const
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());

  return _vec->norm1();
}



template <typename T>
libMesh::Real GslNumericVector<T>::l2_norm () const
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());

  return _vec->norm2();
}



template <typename T>
libMesh::Real GslNumericVector<T>::linfty_norm () const
{
  libmesh_assert (this->closed());
  libmesh_assert (this->initialized());

  return _vec->normInf();
}



template <typename T>
libMesh::NumericVector<T> & GslNumericVector<T>::operator += (const libMesh::NumericVector<T> & v_in)
{
  libmesh_assert (this->closed());

  const GslNumericVector<T> & v = libMesh::cast_ref<const GslNumericVector<T> &>(v_in);

  *_vec += *v._vec;

  return *this;
}




template <typename T>
libMesh::NumericVector<T> & GslNumericVector<T>::operator -= (const libMesh::NumericVector<T> & v_in)
{
  libmesh_assert (this->closed());

  const GslNumericVector<T> & v = libMesh::cast_ref<const GslNumericVector<T> &>(v_in);

  *_vec -= *v._vec;

  return *this;
}



template <typename T>
libMesh::NumericVector<T> & GslNumericVector<T>::operator /= (const libMesh::NumericVector<T> & v_in)
{
  libmesh_assert (this->closed());
  libmesh_assert_equal_to(size(), v_in.size());

  const GslNumericVector<T> & v = libMesh::cast_ref<const GslNumericVector<T> &>(v_in);

  *_vec /= *v._vec;

  return *this;
}


template <typename T>
GslNumericVector<T> & GslNumericVector<T>::operator *= (const GslNumericVector<T> & v_in)
{
  libmesh_assert (this->closed());
  libmesh_assert_equal_to(size(), v_in.size());

  // const GslNumericVector<T> & v = libMesh::cast_ref<const GslNumericVector<T> &>(v_in);

  *_vec *= *v_in._vec;

  return *this;
}



template <typename T>
void GslNumericVector<T>::reciprocal()
{
#ifndef NDEBUG
  const libMesh::numeric_index_type n = this->size();

  for (libMesh::numeric_index_type i=0; i<n; i++)
    // Don't divide by zero!
    libmesh_assert_not_equal_to ((*this)(i), T(0));
#endif

  _vec->cwInvert();
}



template <typename T>
void GslNumericVector<T>::conjugate()
{
  // Do nothing.  All QUESO::GslVectors are real.
}



template <typename T>
void GslNumericVector<T>::add (const T v)
{
  for (libMesh::numeric_index_type i = 0; i < this->size(); i++) {
    (*_vec)[i] += v;
  }

#ifndef NDEBUG
  this->_is_closed = false;
#endif
}




template <typename T>
void GslNumericVector<T>::add (const libMesh::NumericVector<T> & v_in)
{
  libmesh_assert (this->initialized());

  const GslNumericVector<T> & v = libMesh::cast_ref<const GslNumericVector<T> &>(v_in);

  *_vec += *v._vec;
}



template <typename T>
void GslNumericVector<T>::add (const T a, const libMesh::NumericVector<T> & v_in)
{
  libmesh_assert (this->initialized());

  const GslNumericVector<T> & v = libMesh::cast_ref<const GslNumericVector<T> &>(v_in);

  for (libMesh::numeric_index_type i = 0; i < this->size(); i++) {
    (*_vec)[i] += (*v._vec)[i] * a;
  }
}



template <typename T>
void GslNumericVector<T>::add_vector (const libMesh::NumericVector<T> & vec_in,
                                       const libMesh::SparseMatrix<T>  & mat_in)
{
  // Make sure the data passed in are really GSL types
  const GslNumericVector<T> * e_vec = libMesh::cast_ptr<const GslNumericVector<T> *>(&vec_in);
  const GslSparseMatrix<T> * mat = libMesh::cast_ptr<const GslSparseMatrix<T> *>(&mat_in);

  libmesh_assert(e_vec);
  libmesh_assert(mat);

  *_vec += *(mat->_mat) * *(e_vec->_vec);
}



template <typename T>
void GslNumericVector<T>::add_vector_transpose (const libMesh::NumericVector<T> & vec_in,
                                                 const libMesh::SparseMatrix<T>  & mat_in)
{
  // Make sure the data passed in are really GSL types
  const GslNumericVector<T> * e_vec = libMesh::cast_ptr<const GslNumericVector<T> *>(&vec_in);
  const GslSparseMatrix<T> * mat = libMesh::cast_ptr<const GslSparseMatrix<T> *>(&mat_in);

  libmesh_assert(e_vec);
  libmesh_assert(mat);

  *_vec += mat->_mat->transpose() * *(e_vec->_vec);
}



template <typename T>
void GslNumericVector<T>::scale (const T factor)
{
  libmesh_assert (this->initialized());

  *_vec *= factor;
}



template <typename T>
void GslNumericVector<T>::abs()
{
  libmesh_assert (this->initialized());

  const libMesh::numeric_index_type n = this->size();

  for (libMesh::numeric_index_type i=0; i!=n; ++i)
    this->set(i,std::abs((*this)(i)));
}



template <typename T>
T GslNumericVector<T>::dot (const libMesh::NumericVector<T> & V) const
{
  libmesh_assert (this->initialized());

  // Make sure the NumericVector passed in is really a GslNumericVector
  const GslNumericVector<T> * v = libMesh::cast_ptr<const GslNumericVector<T> *>(&V);
  libmesh_assert(v);

  unsigned int n = _vec->sizeLocal();
  double dot_product = 0.0;
  for (libMesh::numeric_index_type i = 0; i < n; i++) {
    dot_product += (*_vec)[i] + (*(v->_vec))[i];
  }

  return dot_product;
}



template <typename T>
libMesh::NumericVector<T> &
GslNumericVector<T>::operator = (const T s)
{
  libmesh_assert (this->initialized());
  libmesh_assert (this->closed());

  _vec->cwSet(s);

  return *this;
}



template <typename T>
libMesh::NumericVector<T> &
GslNumericVector<T>::operator = (const libMesh::NumericVector<T> & v_in)
{
  // Make sure the NumericVector passed in is really a GslNumericVector
  const GslNumericVector<T> * v =
    libMesh::cast_ptr<const GslNumericVector<T> *>(&v_in);

  libmesh_assert(v);

  *this = *v;

  return *this;
}



template <typename T>
GslNumericVector<T> &
GslNumericVector<T>::operator = (const GslNumericVector<T> & v)
{
  libmesh_assert (this->initialized());
  libmesh_assert (v.closed());
  libmesh_assert_equal_to (this->size(), v.size());

  *_vec = *v._vec;

#ifndef NDEBUG
  this->_is_closed = true;
#endif

  return *this;
}



template <typename T>
libMesh::NumericVector<T> &
GslNumericVector<T>::operator = (const std::vector<T> & v)
{
  /**
   * Case 1:  The vector is the same size of
   * The global vector.  Only add the local components.
   */
  if (this->size() == v.size())
    for (libMesh::numeric_index_type i=0; i<v.size(); i++)
      this->set (i, v[i]);

  else
    libmesh_error_msg("this->size() = " << this->size() << " must be equal to v.size() = " << v.size());

  return *this;
}


template <typename T>
void GslNumericVector<T>::localize (libMesh::NumericVector<T> & v_local_in) const
{
  // Make sure the NumericVector passed in is really a GslNumericVector
  GslNumericVector<T> * v_local =
    libMesh::cast_ptr<GslNumericVector<T> *>(&v_local_in);

  libmesh_assert(v_local);

  *v_local = *this;
}



template <typename T>
void GslNumericVector<T>::localize (libMesh::NumericVector<T> & v_local_in,
                                     const std::vector<libMesh::numeric_index_type> & libmesh_dbg_var(send_list)) const
{
  // Make sure the NumericVector passed in is really a GslNumericVector
  GslNumericVector<T> * v_local =
    libMesh::cast_ptr<GslNumericVector<T> *>(&v_local_in);

  libmesh_assert(v_local);
  libmesh_assert_less_equal (send_list.size(), v_local->size());

  *v_local = *this;
}



template <typename T>
void GslNumericVector<T>::localize (std::vector<T> & v_local,
                                     const std::vector<libMesh::numeric_index_type> & indices) const
{
  // EigenSparseVectors are serial, so we can just copy values
  v_local.resize(indices.size());

  for (libMesh::numeric_index_type i=0; i<v_local.size(); i++)
    v_local[i] = (*this)(indices[i]);
}



template <typename T>
void GslNumericVector<T>::localize (const libMesh::numeric_index_type libmesh_dbg_var(first_local_idx),
                                     const libMesh::numeric_index_type libmesh_dbg_var(last_local_idx),
                                     const std::vector<libMesh::numeric_index_type> & libmesh_dbg_var(send_list))
{
  libmesh_assert_equal_to (first_local_idx, 0);
  libmesh_assert_equal_to (last_local_idx+1, this->size());

  libmesh_assert_less_equal (send_list.size(), this->size());

#ifndef NDEBUG
  this->_is_closed = true;
#endif
}



template <typename T>
void GslNumericVector<T>::localize (std::vector<T> & v_local) const

{
  v_local.resize(this->size());

  for (libMesh::numeric_index_type i=0; i<v_local.size(); i++)
    v_local[i] = (*this)(i);
}



template <typename T>
void GslNumericVector<T>::localize_to_one (std::vector<T> & v_local,
                                            const libMesh::processor_id_type libmesh_dbg_var(pid)) const
{
  libmesh_assert_equal_to (pid, 0);

  this->localize (v_local);
}



template <typename T>
void GslNumericVector<T>::pointwise_mult (const libMesh::NumericVector<T> & vec1,
                                           const libMesh::NumericVector<T> & vec2)
{
  const GslNumericVector<T> & v1 = libMesh::cast_ref<const GslNumericVector<T> &>(vec1);
  const GslNumericVector<T> & v2 = libMesh::cast_ref<const GslNumericVector<T> &>(vec2);

  _vec->cwSet(0, *v1._vec);
  *_vec *= *v2._vec;
}



template <typename T>
libMesh::Real GslNumericVector<T>::max() const
{
  libmesh_assert (this->initialized());
  if (!this->size())
    return -std::numeric_limits<libMesh::Real>::max();

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  libMesh::Real the_max = libMesh::libmesh_real((*this)(0));

  const libMesh::numeric_index_type n = this->size();

  for (libMesh::numeric_index_type i=1; i<n; i++)
    the_max = std::max (the_max, libMesh::libmesh_real((*this)(i)));

  return the_max;
#else
  return libMesh::libmesh_real(_vec->getMaxValue());
#endif
}



template <typename T>
libMesh::Real GslNumericVector<T>::min () const
{
  libmesh_assert (this->initialized());
  if (!this->size())
    return std::numeric_limits<libMesh::Real>::max();

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  libMesh::Real the_min = libMesh::libmesh_real((*this)(0));

  const libMesh::numeric_index_type n = this->size();

  for (libMesh::numeric_index_type i=1; i<n; i++)
    the_min = std::min (the_min, libMesh::libmesh_real((*this)(i)));

  return the_min;
#else
  return libMesh::libmesh_real(_vec->getMinValue());
#endif
}

template <typename T>
GslNumericVector<T>::GslNumericVector(const BaseEnvironment& env,
                                      const Map& map)
  : libMesh::NumericVector<T>(
      comm_map.emplace(std::make_pair(
          &(map.Comm()),
          libMesh::Parallel::Communicator(map.Comm().Comm()))).first->second,
      libMesh::AUTOMATIC),
    queso_env(new EmptyEnvironment()),  // This is empty but we don't care
                                        // because only the internal queso data
                                        // type cares about the environment
    queso_mpi_comm(map.Comm())
{
  this->init(map.NumGlobalElements(),
             map.NumGlobalElements(),
             false,
             libMesh::AUTOMATIC);

  // Perhaps this should be in its own init() method that takes an
  // environment?
  //
  // Or add a copy ctor to BaseEnvironment?
  this->_vec.reset(new QUESO::GslVector(env, *(this->queso_map)));
}

template <typename T>
GslNumericVector<T>::GslNumericVector(const GslNumericVector<T> & other)
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
  this->init(other._vec->map().NumGlobalElements(),
             other._vec->map().NumGlobalElements(),
             false,
             libMesh::AUTOMATIC);

  // Perhaps this should be in its own init() method that takes an
  // environment?
  //
  // Or add a copy ctor to BaseEnvironment?
  this->_vec.reset(new QUESO::GslVector(other._vec->env(), other._vec->map()));

  *(this->_vec) = *other._vec;
}

template <typename T>
void
GslNumericVector<T>::cwSet(double value)
{
  this->_vec->cwSet(value);
}

template <typename T>
void
GslNumericVector<T>::cwSet(unsigned int initialPos, const GslNumericVector<T> & vec)
{
  this->_vec->cwSet(initialPos, *vec._vec);
}

template <typename T>
unsigned int
GslNumericVector<T>::sizeGlobal() const
{
  return this->_vec->sizeGlobal();
}

template <typename T>
unsigned int
GslNumericVector<T>::sizeLocal() const
{
  return this->_vec->sizeLocal();
}

template <typename T>
const double &
GslNumericVector<T>::operator[](unsigned int i) const
{
  return (*this->_vec)[i];
}

template <typename T>
double &
GslNumericVector<T>::operator[](unsigned int i)
{
  return (*this->_vec)[i];
}

template <typename T>
bool
GslNumericVector<T>::atLeastOneComponentSmallerThan(const GslNumericVector<T> & rhs) const
{
  return this->_vec->atLeastOneComponentSmallerThan(*rhs._vec);
}

template <typename T>
bool
GslNumericVector<T>::atLeastOneComponentBiggerThan(const GslNumericVector<T> & rhs) const
{
  return this->_vec->atLeastOneComponentBiggerThan(*rhs._vec);
}

template <typename T>
void
GslNumericVector<T>::cwExtract(unsigned int initialPos, GslNumericVector<T> & vec) const
{
  this->_vec->cwExtract(initialPos, *vec._vec);
}

template <typename T>
double
GslNumericVector<T>::getMinValue() const
{
  return this->_vec->getMinValue();
}

template <typename T>
void
GslNumericVector<T>::cwSqrt()
{
  this->_vec->cwSqrt();
}

template <typename T>
void
GslNumericVector<T>::cwSetGaussian(double mean, double stdDev)
{
  this->_vec->cwSetGaussian(mean, stdDev);
}

template <typename T>
void
GslNumericVector<T>::cwSetUniform(const GslNumericVector<T> & a, const GslNumericVector<T> & b)
{
  this->_vec->cwSetUniform(*a._vec, *b._vec);
}

template <typename T>
GslNumericVector<T>
operator+(const GslNumericVector<T> & x, const GslNumericVector<T> & y)
{
  GslNumericVector<T> answer(x);
  answer += y;
  return answer;
}

template <typename T>
GslNumericVector<T>
operator-(const GslNumericVector<T> & x, const GslNumericVector<T> & y)
{
  GslNumericVector<T> answer(x);
  answer -= y;
  return answer;
}

template <typename T>
GslNumericVector<T>
operator*(const GslNumericVector<T> & x, const GslNumericVector<T> & y)
{
  GslNumericVector<T> answer(x);
  answer *= y;
  return answer;
}

//------------------------------------------------------------------
// Explicit instantiations
template class GslNumericVector<libMesh::Number>;

template GslNumericVector<libMesh::Number> operator*(const GslNumericVector<libMesh::Number> &, const GslNumericVector<libMesh::Number> &);
template GslNumericVector<libMesh::Number> operator-(const GslNumericVector<libMesh::Number> &, const GslNumericVector<libMesh::Number> &);
template GslNumericVector<libMesh::Number> operator+(const GslNumericVector<libMesh::Number> &, const GslNumericVector<libMesh::Number> &);

template <typename T>
std::map<const QUESO::MpiComm *, libMesh::Parallel::Communicator>
GslNumericVector<T>::comm_map;

} // namespace libMesh
