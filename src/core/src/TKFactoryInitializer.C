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

#include <queso/TKFactoryInitializer.h>

#include <queso/TKFactoryRandomWalk.h>
#include <queso/TKFactoryMALA.h>
#include <queso/TKFactoryLogitRandomWalk.h>
#include <queso/TKFactoryStochasticNewton.h>
#include <queso/ScaledCovMatrixTKGroup.h>
#include <queso/TransformedScaledCovMatrixTKGroup.h>
#include <queso/MetropolisAdjustedLangevinTK.h>
#include <queso/HessianCovMatricesTKGroup.h>

#include <queso/GslVector.h>
#include <queso/GslNumericVector.h>
#include <queso/GslMatrix.h>
#include <queso/GslSparseMatrix.h>

namespace QUESO
{

template <typename V, typename M>
TKFactoryInitializer<V, M>::TKFactoryInitializer()
{
  // Instantiate all the transition kernel factories
  static TKFactoryRandomWalk<ScaledCovMatrixTKGroup<V, M>, V, M> tk_factory_random_walk("random_walk");
  static TKFactoryLogitRandomWalk<TransformedScaledCovMatrixTKGroup<V, M>, V, M> tk_factory_logit_random_walk("logit_random_walk");
  static TKFactoryStochasticNewton<HessianCovMatricesTKGroup<V, M>, V, M> tk_factory_stochastic_newton("stochastic_newton");
  static TKFactoryMALA<MetropolisAdjustedLangevinTK<V, M>, V, M> tk_factory_mala("mala");
}

template <typename V, typename M>
TKFactoryInitializer<V, M>::~TKFactoryInitializer()
{
  // Do nothing
}

template class TKFactoryInitializer<GslVector, GslMatrix>;
template class TKFactoryInitializer<GslNumericVector<libMesh::Number>, GslSparseMatrix<libMesh::Number> >;

} // namespace QUESO
