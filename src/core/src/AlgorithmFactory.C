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

#include <queso/GslVector.h>
#include <queso/GslNumericVector.h>
#include <queso/GslMatrix.h>
#include <queso/GslSparseMatrix.h>
#include <libmesh/eigen_sparse_vector.h>
#include <libmesh/eigen_sparse_matrix.h>
#include <queso/AlgorithmFactory.h>

namespace QUESO
{

template <>
std::map<std::string, Factory<Algorithm<GslVector, GslMatrix> > *> &
Factory<Algorithm<GslVector, GslMatrix> >::factory_map()
{
  static std::map<std::string, Factory<Algorithm<GslVector, GslMatrix> > *> _factory_map;

  return _factory_map;
}

template <>
std::map<std::string, Factory<Algorithm<GslNumericVector<libMesh::Number>, GslSparseMatrix<libMesh::Number> > > *> &
Factory<Algorithm<GslNumericVector<libMesh::Number>, GslSparseMatrix<libMesh::Number> > >::factory_map()
{
  static std::map<std::string, Factory<Algorithm<GslNumericVector<libMesh::Number>, GslSparseMatrix<libMesh::Number> > > *> _factory_map;

  return _factory_map;
}

template <>
std::map<std::string, Factory<Algorithm<libMesh::EigenSparseVector<libMesh::Number>, libMesh::EigenSparseMatrix<libMesh::Number> > > *> &
Factory<Algorithm<libMesh::EigenSparseVector<libMesh::Number>, libMesh::EigenSparseMatrix<libMesh::Number> > >::factory_map()
{
  static std::map<std::string, Factory<Algorithm<libMesh::EigenSparseVector<libMesh::Number>, libMesh::EigenSparseMatrix<libMesh::Number> > > *> _factory_map;

  return _factory_map;
}

template <>
const BaseEnvironment * AlgorithmFactory<GslVector, GslMatrix>::m_env = NULL;

template <>
const BaseEnvironment * AlgorithmFactory<GslNumericVector<libMesh::Number>, GslSparseMatrix<libMesh::Number> >::m_env = NULL;

template <>
const BaseEnvironment * AlgorithmFactory<libMesh::EigenSparseVector<libMesh::Number>, libMesh::EigenSparseMatrix<libMesh::Number> >::m_env = NULL;

template <>
const BaseTKGroup<GslVector, GslMatrix> * AlgorithmFactory<GslVector, GslMatrix>::m_tk = NULL;

template <>
const BaseTKGroup<GslNumericVector<libMesh::Number>, GslSparseMatrix<libMesh::Number> > * AlgorithmFactory<GslNumericVector<libMesh::Number>, GslSparseMatrix<libMesh::Number> >::m_tk = NULL;

template <>
const BaseTKGroup<libMesh::EigenSparseVector<libMesh::Number>, libMesh::EigenSparseMatrix<libMesh::Number> > * AlgorithmFactory<libMesh::EigenSparseVector<libMesh::Number>, libMesh::EigenSparseMatrix<libMesh::Number> >::m_tk = NULL;

} // namespace QUESO
