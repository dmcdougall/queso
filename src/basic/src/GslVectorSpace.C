//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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

#include <queso/VectorSpace.h>
#include <queso/GslMatrix.h>
#include <queso/GslNumericVector.h>
#include <queso/GslSparseMatrix.h>
#include <libmesh/eigen_sparse_vector.h>
#include <libmesh/eigen_sparse_matrix.h>

namespace QUESO {

template <>
Map*
VectorSpace<GslVector, GslMatrix>::newMap()
{
  return new Map(m_dimGlobal,0,m_env.selfComm());
}

template <>
Map*
VectorSpace<GslNumericVector<libMesh::Number>, GslSparseMatrix<libMesh::Number> >::newMap()
{
  return new Map(m_dimGlobal,0,m_env.selfComm());
}

template <>
Map*
VectorSpace<libMesh::EigenSparseVector<libMesh::Number>, libMesh::EigenSparseMatrix<libMesh::Number> >::newMap()
{
  return new Map(m_dimGlobal,0,m_env.selfComm());
}


template<>
GslVector*
VectorSpace<GslVector,GslMatrix>::newVector() const
{
  return new GslVector(m_env,*m_map);
}

template<>
GslNumericVector<libMesh::Number>*
VectorSpace<GslNumericVector<libMesh::Number>,GslSparseMatrix<libMesh::Number> >::newVector() const
{
  return new GslNumericVector<libMesh::Number>(m_env,*m_map);
}

template<>
libMesh::EigenSparseVector<libMesh::Number>*
VectorSpace<libMesh::EigenSparseVector<libMesh::Number>,libMesh::EigenSparseMatrix<libMesh::Number> >::newVector() const
{
  return new libMesh::EigenSparseVector<libMesh::Number>(m_env,*m_map);
}

template<>
GslVector*
VectorSpace<GslVector,GslMatrix>::newVector(double value) const
{
  return new GslVector(m_env,*m_map,value);
}

template<>
GslNumericVector<libMesh::Number>*
VectorSpace<GslNumericVector<libMesh::Number>,GslSparseMatrix<libMesh::Number> >::newVector(double value) const
{
  GslNumericVector<libMesh::Number> * result =
    new GslNumericVector<libMesh::Number>(m_env,*m_map);

  result->cwSet(value);

  return result;
}

template<>
libMesh::EigenSparseVector<libMesh::Number>*
VectorSpace<libMesh::EigenSparseVector<libMesh::Number>,libMesh::EigenSparseMatrix<libMesh::Number> >::newVector(double value) const
{
  libMesh::EigenSparseVector<libMesh::Number> * result =
    new libMesh::EigenSparseVector<libMesh::Number>(m_env,*m_map);

  result->cwSet(value);

  return result;
}

template<>
libMesh::EigenSparseVector<libMesh::Number>*
VectorSpace<libMesh::EigenSparseVector<libMesh::Number>,libMesh::EigenSparseVector<libMesh::Number> >::newVector(double value) const
{
  libMesh::EigenSparseVector<libMesh::Number> * result =
    new libMesh::EigenSparseVector<libMesh::Number>(m_env,*m_map);

  result->cwSet(value);

  return result;
}

template<>
GslMatrix*
VectorSpace<GslVector,GslMatrix>::newMatrix() const
{
  return new GslMatrix(m_env,*m_map,this->dimGlobal());
}

template<>
GslSparseMatrix<libMesh::Number>*
VectorSpace<GslNumericVector<libMesh::Number>, GslSparseMatrix<libMesh::Number> >::newMatrix() const
{
  return new GslSparseMatrix<libMesh::Number>(m_env,*m_map,this->dimGlobal());
}

template<>
libMesh::EigenSparseMatrix<libMesh::Number>*
VectorSpace<libMesh::EigenSparseVector<libMesh::Number>, libMesh::EigenSparseMatrix<libMesh::Number> >::newMatrix() const
{
  return new libMesh::EigenSparseMatrix<libMesh::Number>(m_env,*m_map,this->dimGlobal());
}

template<>
GslMatrix*
VectorSpace<GslVector,GslMatrix>::newDiagMatrix(double diagValue) const
{
  return new GslMatrix(m_env,*m_map,diagValue);
}

template<>
GslSparseMatrix<libMesh::Number> *
VectorSpace<GslNumericVector<libMesh::Number>,GslSparseMatrix<libMesh::Number> >::newDiagMatrix(double diagValue) const
{
  return new GslSparseMatrix<libMesh::Number>(m_env,*m_map,diagValue);
}

}  // End namespace QUESO
