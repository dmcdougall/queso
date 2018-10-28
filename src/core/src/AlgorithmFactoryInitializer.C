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

#include <queso/AlgorithmFactory.h>
#include <queso/AlgorithmFactoryInitializer.h>
#include <queso/GslVector.h>
#include <queso/GslNumericVector.h>
#include <queso/GslMatrix.h>
#include <queso/GslSparseMatrix.h>
#include <libmesh/eigen_sparse_vector.h>
#include <libmesh/eigen_sparse_matrix.h>

namespace QUESO
{

template <typename V, typename M>
AlgorithmFactoryInitializer<V, M>::AlgorithmFactoryInitializer()
{
  // Instantiate all the algorithm factories
  static AlgorithmFactoryImp<Algorithm, V, M> random_walk_alg("random_walk");
  static AlgorithmFactoryImp<Algorithm, V, M> logit_random_walk_alg("logit_random_walk");

}

template <typename V, typename M>
AlgorithmFactoryInitializer<V, M>::~AlgorithmFactoryInitializer()
{
  // Do nothing
}

template class AlgorithmFactoryInitializer<GslVector, GslMatrix>;
template class AlgorithmFactoryInitializer<GslNumericVector<libMesh::Number>, GslSparseMatrix<libMesh::Number> >;
template class AlgorithmFactoryInitializer<libMesh::EigenSparseVector<libMesh::Number>, libMesh::EigenSparseMatrix<libMesh::Number> >;

} // namespace QUESO
