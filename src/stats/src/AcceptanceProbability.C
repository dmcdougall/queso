//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
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
#include <queso/GslMatrix.h>
#include <queso/AcceptanceProbability.h>

namespace QUESO {

template <class V, class M>
AcceptanceProbability<V, M>::AcceptanceProbability()
{
}

template <class V, class M>
AcceptanceProbability<V, M>::~AcceptanceProbability()
{
}

template <class V, class M>
double
AcceptanceProbability<V, M>::evaluate(const MarkovChainPositionData<V, M> & x,
                                      const MarkovChainPositionData<V, M> & y,
                                      unsigned int xStageId,
                                      unsigned int yStageId,
                                      double & alphaQuotient)
{
  return 0.0;
}

template class AcceptanceProbability<GslVector, GslMatrix>;
}  // End namespace QUESO
