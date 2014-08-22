//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#include <cmath>
#include <gsl/gsl_sort_vector.h>
#include <queso/GslVectorImplementation.h>

namespace QUESO {

GslVectorImplementation::GslVectorImplementation(const Map& map)
  : m_vec(gsl_vector_calloc(map.NumGlobalElements()))
{
  if (m_vec == NULL) {
    std::cerr << "GslVectorImplementation::constructor(1): null vector "
              << "generated" << std::endl;
  }

  if (m_vec->size != (unsigned int) map.NumMyElements()) {
    std::cerr << "GslVectorImplementation::constructor(1): incompatible local "
              << "vec size" << std::endl;
  }

  if (m_vec->size != (unsigned int) map.NumGlobalElements()) {
    std::cerr << "GslVectorImplementation::constructor(1): incompatible "
              << "global vec size" << std::endl;
  }

  if (m_vec->size != (unsigned int) map.NumMyElements()) {
    std::cerr << "GslVectorImplementation::constructor(1): incompatible own "
              << "vec size" << std::endl;
  }
}

GslVectorImplementation::GslVectorImplementation(const Map& map, double value)
  : m_vec(gsl_vector_calloc(map.NumGlobalElements()))
{
  if (m_vec == NULL) {
    std::cerr << "GslVectorImplementation::constructor(2): null vector "
              << "generated" << std::endl;
  }

  if (m_vec->size != (unsigned int) map.NumMyElements()) {
    std::cerr << "GslVectorImplementation::constructor(2): incompatible "
              << "local vec size" << std::endl;
  }

  if (m_vec->size != (unsigned int) map.NumGlobalElements()) {
    std::cerr << "GslVectorImplementation::constructor(2): incompatible "
              << "global vec size" << std::endl;
  }

  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    (*this)[i] = value;
  }

  if (m_vec->size != (unsigned int) map.NumMyElements()) {
    std::cerr << "GslVectorImplementation::constructor(2): incompatible own "
              << "vec size" << std::endl;
  }
}

GslVectorImplementation::GslVectorImplementation(double d1, double d2,
    const Map& map)
  : m_vec(gsl_vector_calloc(map.NumGlobalElements()))
{
  if (m_vec == NULL) {
    std::cerr << "GslVectorImplementation::constructor(3), linspace: null "
              << "vector generated" << std::endl;
  }

  if (m_vec->size != (unsigned int) map.NumMyElements()) {
    std::cerr << "GslVectorImplementation::constructor(3): incompatible "
              << "local vec size" << std::endl;
  }

  if (m_vec->size != (unsigned int) map.NumGlobalElements()) {
    std::cerr << "GslVectorImplementation::constructor(3): incompatible "
              << "global vec size" << std::endl;
  }

  for (unsigned int i = 0; i < m_vec->size; ++i) {
    double alpha = (double) i / ((double) m_vec->size - 1.0);
    (*this)[i] = (1.0 - alpha) * d1 + alpha * d2;
  }

  if (m_vec->size != (unsigned int) map.NumMyElements()) {
    std::cerr << "GslVectorImplementation::constructor(3): incompatible "
              << "own vec size" << std::endl;
  }
}

GslVectorImplementation::GslVectorImplementation(
    const GslVectorImplementation& v, double start, double end)
  : m_vec(gsl_vector_calloc(v.sizeLocal()))
{
  if (m_vec == NULL) {
    std::cerr << "GslVectorImplementation::constructor(4), linspace: null "
              << "vector generated" << std::endl;
  }

  for (unsigned int i = 0; i < m_vec->size; ++i) {
    double alpha = (double) i / ((double) m_vec->size - 1.0);
    (*this)[i] = (1.0 - alpha) * start + alpha * end;
  }
}

GslVectorImplementation::GslVectorImplementation(
    const GslVectorImplementation& v)
  : m_vec(gsl_vector_calloc(v.sizeLocal()))
{
  // FIX THIS
  if (m_vec == NULL) {
    std::cerr << "GslVectorImplementation::constructor(5), copy: null vector "
              << "generated" << std::endl;
  }

  this->copy(v);
}

GslVectorImplementation::~GslVectorImplementation()
{
  if (m_vec) {
    gsl_vector_free(m_vec);
  }
}

GslVectorImplementation&
GslVectorImplementation::operator=(const GslVectorImplementation& rhs)
{
  // FIX THIS
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();
  if (size1 != size2) {
    std::cerr << "GslVectorImplementation::operator=(): sizes are not "
              << "compatible" << std::endl;
  }
  this->copy(rhs);
  return *this;
}

GslVectorImplementation&
GslVectorImplementation::operator*=(double a)
{
  int iRC;
  iRC = gsl_vector_scale(m_vec, a);
  if (iRC) {
    std::cerr << "GslVectorImplementation::operator*=(): failed" << std::endl;
  }
  return *this;
}

GslVectorImplementation&
GslVectorImplementation::operator+=(const GslVectorImplementation& rhs)
{
  int iRC;
  iRC = gsl_vector_add(m_vec, rhs.m_vec);
  if (iRC) {
    std::cerr << "GslVectorImplementation::operator+=(): failed" << std::endl;
  }
  return *this;
}

GslVectorImplementation&
GslVectorImplementation::operator-=(const GslVectorImplementation& rhs)
{
  int iRC;
  iRC = gsl_vector_sub(m_vec, rhs.m_vec);
  if (iRC) {
    std::cerr << "GslVectorImplementation::operator-=(): failed" << std::endl;
  }

  return *this;
}

double&
GslVectorImplementation::operator[](unsigned int i)
{
  return *gsl_vector_ptr(m_vec,i);
}

const double&
GslVectorImplementation::operator[](unsigned int i) const
{
  return *gsl_vector_const_ptr(m_vec,i);
}

void
GslVectorImplementation::copy(const GslVectorImplementation& src)
{
  int iRC;
  iRC = gsl_vector_memcpy(this->m_vec, src.m_vec);
  if (iRC) {
    std::cerr << "GslVectorImplementation::copy(): failed" << std::endl;
  }

  return;
}

unsigned int
GslVectorImplementation::sizeLocal() const
{
  return m_vec->size;
}

unsigned int
GslVectorImplementation::sizeGlobal() const
{
  return m_vec->size;
}

void
GslVectorImplementation::sort()
{
  gsl_sort_vector(m_vec);
}

gsl_vector*
GslVectorImplementation::data() const
{
  return m_vec;
}

double
GslVectorImplementation::getMaxValue() const
{
  return gsl_vector_max(m_vec);
}

double
GslVectorImplementation::getMinValue() const
{
  return gsl_vector_min(m_vec);
}

int
GslVectorImplementation::getMaxValueIndex() const
{
  return gsl_vector_max_index(m_vec);
}

int
GslVectorImplementation::getMinValueIndex() const
{
  return gsl_vector_min_index(m_vec);
}

}  // End namespace QUESO
