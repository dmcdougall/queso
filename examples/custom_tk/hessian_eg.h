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

#ifndef QUESO_CUSTOM_TK_HESS_H
#define QUESO_CUSTOM_TK_HESS_H

#include <queso/TKGroup.h>
#include <queso/VectorRV.h>
#include <queso/ScalarFunctionSynchronizer.h>

namespace QUESO {

class GslVector;
class GslMatrix;

template <class V = GslVector, class M = GslMatrix>
class MyTransitionKernel : public BaseTKGroup<V,M> {
public:
  MyTransitionKernel(const char * prefix,
                     const VectorSpace<V, M> & vectorSpace,
                     const std::vector<double> & scales,
                     const M & covMatrix);

  ~MyTransitionKernel();

  bool symmetric() const;
  const GaussianVectorRV<V, M> & rv(unsigned int stageId) const;
  const GaussianVectorRV<V, M> & rv(const std::vector<unsigned int> & stageIds);
  virtual const GaussianVectorRV<V, M> & rv(const V & position) const;
  void updateLawCovMatrix(const M & covMatrix);
  bool setPreComputingPosition(const V & position, unsigned int stageId);
  void clearPreComputingPositions();
  virtual unsigned int set_dr_stage(unsigned int stageId);
  virtual void update_tk();
  void print(std::ostream & os) const;
private:
  void setRVsWithZeroMean();
  using BaseTKGroup<V,M>::m_env;
  using BaseTKGroup<V,M>::m_prefix;
  using BaseTKGroup<V,M>::m_vectorSpace;
  using BaseTKGroup<V,M>::m_scales;
  using BaseTKGroup<V,M>::m_preComputingPositions;
  using BaseTKGroup<V,M>::m_rvs;

  M m_originalCovMatrix;
};

}  // End namespace QUESO

#endif // QUESO_CUSTOM_TK_HESS_H
