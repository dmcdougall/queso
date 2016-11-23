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

#include "hessian_eg.h"

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/GaussianJointPdf.h>
#include <queso/TKFactoryRandomWalk.h>

namespace QUESO {

template<class V, class M>
MyTransitionKernel<V,M>::MyTransitionKernel(
  const char*                    prefix,
  const VectorSpace<V,M>& vectorSpace, // FIX ME: vectorSubset ???
  const std::vector<double>&     scales,
  const M&                       covMatrix)
  :
  BaseTKGroup<V,M>(prefix,vectorSpace,scales),
  m_originalCovMatrix    (covMatrix)
{
  std::cout << "Hello from MyTransitionKernel!" << std::endl;
  setRVsWithZeroMean();
}

template<class V, class M>
MyTransitionKernel<V,M>::~MyTransitionKernel()
{
}

template<class V, class M>
bool
MyTransitionKernel<V,M>::symmetric() const
{
  return false;
}

template<class V, class M>
const GaussianVectorRV<V,M>&
MyTransitionKernel<V,M>::rv(unsigned int stageId) const
{
  queso_require_not_equal_to_msg(m_rvs.size(), 0, "m_rvs.size() = 0");
  queso_require_msg(m_rvs[0], "m_rvs[0] == NULL");
  queso_require_greater_msg(m_preComputingPositions.size(), stageId, "m_preComputingPositions.size() <= stageId");
  queso_require_msg(m_preComputingPositions[stageId], "m_preComputingPositions[stageId] == NULL");

  if ((m_env.subDisplayFile()        ) &&
      (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "In MyTransitionKernel<V,M>::rv1()"
                            << ", stageId = " << stageId
                            << ": about to call m_rvs[0]->updateLawExpVector()"
                            << ", vector = " << *m_preComputingPositions[stageId] // FIX ME: might demand parallelism
                            << std::endl;
  }

  GaussianVectorRV<V, M> * gaussian_rv = dynamic_cast<GaussianVectorRV<V, M> * >(m_rvs[0]);

  gaussian_rv->updateLawExpVector(*m_preComputingPositions[stageId]);

  return (*gaussian_rv);
}

template<class V, class M>
const GaussianVectorRV<V,M>&
MyTransitionKernel<V,M>::rv(const std::vector<unsigned int>& stageIds)
{
  queso_require_greater_equal_msg(m_rvs.size(), stageIds.size(), "m_rvs.size() < stageIds.size()");
  queso_require_msg(m_rvs[stageIds.size()-1], "m_rvs[stageIds.size()-1] == NULL");
  queso_require_greater_msg(m_preComputingPositions.size(), stageIds[0], "m_preComputingPositions.size() <= stageIds[0]");
  queso_require_msg(m_preComputingPositions[stageIds[0]], "m_preComputingPositions[stageIds[0]] == NULL");

  if ((m_env.subDisplayFile()        ) &&
      (m_env.displayVerbosity() >= 10)) {
    *m_env.subDisplayFile() << "In MyTransitionKernel<V,M>::rv2()"
                            << ", stageIds.size() = " << stageIds.size()
                            << ", stageIds[0] = "     << stageIds[0]
                            << ": about to call m_rvs[stageIds.size()-1]->updateLawExpVector()"
                            << ", vector = " << *m_preComputingPositions[stageIds[0]] // FIX ME: might demand parallelism
                            << std::endl;
  }

  GaussianVectorRV<V, M> * gaussian_rv = dynamic_cast<GaussianVectorRV<V, M> * >(m_rvs[stageIds.size()-1]);

  gaussian_rv->updateLawExpVector(*m_preComputingPositions[stageIds[0]]);

  return (*gaussian_rv);
}

template <class V, class M>
const GaussianVectorRV<V, M> &
MyTransitionKernel<V, M>::rv(const V & position) const
{
  queso_require_not_equal_to_msg(m_rvs.size(), 0, "m_rvs.size() = 0");
  queso_require_msg(m_rvs[0], "m_rvs[0] == NULL");
  //queso_require_greater_msg(m_preComputingPositions.size(), this->m_stageId, "m_preComputingPositions.size() <= stageId");
  //queso_require_msg(m_preComputingPositions[this->m_stageId], "m_preComputingPositions[stageId] == NULL");

  // QUESO relies on the user not touching the RV's covariance matrix?
  GaussianVectorRV<V, M> * gaussian_rv = dynamic_cast<GaussianVectorRV<V, M> * >(m_rvs[this->m_stageId]);

  gaussian_rv->updateLawExpVector(position);

  return (*gaussian_rv);
}

template<class V, class M>
void
MyTransitionKernel<V,M>::updateLawCovMatrix(const M& covMatrix)
{
  for (unsigned int i = 0; i < m_scales.size(); ++i) {
    double factor = 1./m_scales[i]/m_scales[i];
    if ((m_env.subDisplayFile()        ) &&
        (m_env.displayVerbosity() >= 10)) {
      *m_env.subDisplayFile() << "In MyTransitionKernel<V,M>::updateLawCovMatrix()"
                              << ", m_scales.size() = " << m_scales.size()
                              << ", i = "               << i
                              << ", m_scales[i] = "     << m_scales[i]
                              << ", factor = "          << factor
                              << ": about to call m_rvs[i]->updateLawCovMatrix()"
                              << ", covMatrix = \n" << factor*covMatrix // FIX ME: might demand parallelism
                              << std::endl;
    }
    dynamic_cast<GaussianVectorRV<V, M> * >(m_rvs[i])->updateLawCovMatrix(factor*covMatrix);
  }
}

template<class V, class M>
bool
MyTransitionKernel<V,M>::setPreComputingPosition(const V& position, unsigned int stageId)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering MyTransitionKernel<V,M>::setPreComputingPosition()"
                           << ": position = " << position
                           << ", stageId = "  << stageId
                           << std::endl;
  }

  BaseTKGroup<V,M>::setPreComputingPosition(position,stageId);
  //setRVsWithZeroMean();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In MyTransitionKernel<V,M>::setPreComputingPosition()"
                           << ", position = "        << position
                           << ", stageId = "         << stageId
                           << ": preComputingPos = " << *m_preComputingPositions[stageId];
    if (stageId < m_scales.size()) {
      *m_env.subDisplayFile() << ", factor = " << 1./m_scales[stageId]/m_scales[stageId];
    }
    if (stageId < m_rvs.size()) {
      const GaussianJointPdf<V,M>* pdfPtr = dynamic_cast< const GaussianJointPdf<V,M>* >(&(m_rvs[stageId]->pdf()));
      *m_env.subDisplayFile() << ", rvCov = " << pdfPtr->lawCovMatrix(); // FIX ME: might demand parallelism
    }
    *m_env.subDisplayFile() << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving MyTransitionKernel<V,M>::setPreComputingPosition()"
                           << ": position = " << position
                           << ", stageId = "  << stageId
                           << std::endl;
  }

  return true;
}

template<class V, class M>
void
MyTransitionKernel<V,M>::clearPreComputingPositions()
{
  BaseTKGroup<V,M>::clearPreComputingPositions();
}

template <class V, class M>
unsigned int
MyTransitionKernel<V, M>::set_dr_stage(unsigned int stageId)
{
  unsigned int old_stageId = this->m_stageId;
  this->m_stageId = stageId;
  return old_stageId;
}

template<class V, class M>
void
MyTransitionKernel<V,M>::setRVsWithZeroMean()
{
  queso_require_not_equal_to_msg(m_rvs.size(), 0, "m_rvs.size() = 0");

  queso_require_equal_to_msg(m_rvs.size(), m_scales.size(), "m_rvs.size() != m_scales.size()");

  for (unsigned int i = 0; i < m_scales.size(); ++i) {
    double factor = 1./m_scales[i]/m_scales[i];
    queso_require_msg(!(m_rvs[i]), "m_rvs[i] != NULL");
    m_rvs[i] = new GaussianVectorRV<V,M>(m_prefix.c_str(),
                                         *m_vectorSpace,
                                         m_vectorSpace->zeroVector(),
                                         factor*m_originalCovMatrix);
  }
}

template <class V, class M>
void
MyTransitionKernel<V, M>::update_tk()
{
  std::cout << "QUESO called `update_tk'" << std::endl;

  GslMatrix new_matrix(m_vectorSpace->zeroVector());
  new_matrix(0, 0) = 1.0;

  // Modify transition kernel's random variable here
  for (unsigned int i = 0; i < m_scales.size(); ++i) {
    delete m_rvs[i];
    double factor = 1.0 / m_scales[i] / m_scales[i];
    m_rvs[i] = new GaussianVectorRV<V, M>(m_prefix.c_str(),
        *m_vectorSpace,
        m_vectorSpace->zeroVector(),
        factor*new_matrix);
    std::cout << "Updating rv for dr stage " << i << " with cov matrix: " << std::endl;
    std::cout << factor * new_matrix << std::endl;
  }
}

template<class V, class M>
void
MyTransitionKernel<V,M>::print(std::ostream& os) const
{
  BaseTKGroup<V,M>::print(os);
}

// Explicit instantiation of the template
template class MyTransitionKernel<GslVector, GslMatrix>;

// Register this TK with the appropriate factory
TKFactoryRandomWalk<MyTransitionKernel<GslVector, GslMatrix> > tk_factory_mytk("my_tk");

}  // End namespace QUESO
