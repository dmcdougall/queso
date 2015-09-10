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
#include <queso/VectorRV.h>
#include <queso/TKGroup.h>
#include <queso/ScalarFunctionSynchronizer.h>
#include <queso/MetropolisHastingsSGOptions.h>
#include <queso/MarkovChainPositionData.h>

#include <queso/SequenceGenerator.h>

namespace QUESO {

template<class V, class M>
SequenceGenerator<V,M>::SequenceGenerator(const char * prefix,
                                          const MhOptionsValues * alternativeOptionsValues,
                                          const BaseVectorRV<V, M> & sourceRv,
                                          const V & initialPosition,
                                          const M * inputProposalCovMatrix)
: m_env(sourceRv.env()),
  m_vectorSpace(sourceRv.imageSet().vectorSpace()),
  m_targetPdf(sourceRv.pdf()),
  m_initialPosition(initialPosition),
  m_initialProposalCovMatrix(m_vectorSpace.zeroVector()),
  m_nullInputProposalCovMatrix(inputProposalCovMatrix == NULL),
  m_targetPdfSynchronizer(new ScalarFunctionSynchronizer<V,M>(m_targetPdf,m_initialPosition)),
  m_tk(NULL),
  m_idsOfUniquePositions(0),
  m_logTargets(0),
  m_alphaQuotients(0),
  m_lastChainSize(0),
  m_lastMean(NULL),
  m_lastAdaptedCovMatrix(NULL),
  m_optionsObj(alternativeOptionsValues),
  m_computeInitialPriorAndLikelihoodValues(true),
  m_initialLogPriorValue(0.),
  m_initialLogLikelihoodValue(0.),
  m_userDidNotProvideOptions(false)
{
  if (inputProposalCovMatrix != NULL) {
    m_initialProposalCovMatrix = *inputProposalCovMatrix;
  }

  queso_require_equal_to_msg(sourceRv.imageSet().vectorSpace().dimLocal(), initialPosition.sizeLocal(), "'sourceRv' and 'initialPosition' should have equal dimensions");

  if (inputProposalCovMatrix) {
    queso_require_equal_to_msg(sourceRv.imageSet().vectorSpace().dimLocal(), inputProposalCovMatrix->numRowsLocal(), "'sourceRv' and 'inputProposalCovMatrix' should have equal dimensions");
    queso_require_equal_to_msg(inputProposalCovMatrix->numCols(), inputProposalCovMatrix->numRowsGlobal(), "'inputProposalCovMatrix' should be a square matrix");
  }
}

template<class V, class M>
SequenceGenerator<V, M>::~SequenceGenerator()
{
  if (m_lastAdaptedCovMatrix) {
    delete m_lastAdaptedCovMatrix;
  }
  if (m_lastMean) {
    delete m_lastMean;
  }
  m_lastChainSize = 0;
  m_alphaQuotients.clear();
  m_logTargets.clear();
  m_idsOfUniquePositions.clear();

  if (m_tk) {
    delete m_tk;
  }
  if (m_targetPdfSynchronizer) {
    delete m_targetPdfSynchronizer;
  }

  // Only delete if the user didn't provide the options
  // I.e., if the user created their options object, then they are resonsible
  // for freeing it.
  if (m_optionsObj && m_userDidNotProvideOptions) {
    delete m_optionsObj;
  }
}

template<class V, class M>
double
SequenceGenerator<V,M>::alpha(const MarkovChainPositionData<V> & x,
                              const MarkovChainPositionData<V> & y,
                              unsigned int xStageId,
                              unsigned int yStageId,
                              double * alphaQuotientPtr)
{
  double alphaQuotient = 0.;
  if ((x.outOfTargetSupport() == false) &&
      (y.outOfTargetSupport() == false)) {
    if ((x.logTarget() == -INFINITY) ||
        (x.logTarget() ==  INFINITY) ||
        ( (boost::math::isnan)(x.logTarget())      )) {
      std::cerr << "WARNING In SequenceGenerator<V,M>::alpha(x,y)"
                << ", worldRank "       << m_env.worldRank()
                << ", fullRank "        << m_env.fullRank()
                << ", subEnvironment "  << m_env.subId()
                << ", subRank "         << m_env.subRank()
                << ", inter0Rank "      << m_env.inter0Rank()
                << ": x.logTarget() = " << x.logTarget()
                << ", x.values() = "    << x.vecValues()
                << ", y.values() = "    << y.vecValues()
                << std::endl;
    }
    else if ((y.logTarget() == -INFINITY           ) ||
             (y.logTarget() ==  INFINITY           ) ||
             ( (boost::math::isnan)(y.logTarget()) )) {
      std::cerr << "WARNING In SequenceGenerator<V,M>::alpha(x,y)"
                << ", worldRank "       << m_env.worldRank()
                << ", fullRank "        << m_env.fullRank()
                << ", subEnvironment "  << m_env.subId()
                << ", subRank "         << m_env.subRank()
                << ", inter0Rank "      << m_env.inter0Rank()
                << ": y.logTarget() = " << y.logTarget()
                << ", x.values() = "    << x.vecValues()
                << ", y.values() = "    << y.vecValues()
                << std::endl;
    }
    else {
      double yLogTargetToUse = y.logTarget();

      if (m_tk->symmetric()) {
        alphaQuotient = std::exp(yLogTargetToUse - x.logTarget());

        if ((m_env.subDisplayFile()                   ) &&
            (m_env.displayVerbosity() >= 3            ) &&
            (m_optionsObj->m_totallyMute == false)) {
          *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::alpha(x,y)"
                                 << ": symmetric proposal case"
                                 << ", x = "               << x.vecValues()
                                 << ", y = "               << y.vecValues()
                                 << ", yLogTargetToUse = " << yLogTargetToUse
                                 << ", x.logTarget() = "   << x.logTarget()
                                 << ", alpha = "           << alphaQuotient
                                 << std::endl;
        }
      }
      else {
        double qyx = m_tk->rv(yStageId).pdf().lnValue(x.vecValues(),NULL,NULL,NULL,NULL);
        double qxy = m_tk->rv(xStageId).pdf().lnValue(y.vecValues(),NULL,NULL,NULL,NULL);
        alphaQuotient = std::exp(yLogTargetToUse +
                                 qyx -
                                 x.logTarget() -
                                 qxy);
        if ((m_env.subDisplayFile()                   ) &&
            (m_env.displayVerbosity() >= 3            ) &&
            (m_optionsObj->m_totallyMute == false)) {
          *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::alpha(x,y)"
                                 << ": asymmetric proposal case"
                                 << ", xStageId = "        << xStageId
                                 << ", yStageId = "        << yStageId
                                 << ", x = "               << x.vecValues()
                                 << ", y = "               << y.vecValues()
                                 << ", yLogTargetToUse = " << yLogTargetToUse
                                 << ", q(y,x) = "          << qyx
                                 << ", x.logTarget() = "   << x.logTarget()
                                 << ", q(x,y) = "          << qxy
                                 << ", alpha = "           << alphaQuotient
                                 << std::endl;
        }
      }
    } // protection logic against logTarget values
  }
  else {
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 10           ) &&
        (m_optionsObj->m_totallyMute == false)) {
      *m_env.subDisplayFile() << "In SequenceGenerator<V,M>::alpha(x,y)"
                             << ": x.outOfTargetSupport = " << x.outOfTargetSupport()
                             << ", y.outOfTargetSupport = " << y.outOfTargetSupport()
                             << std::endl;
    }
  }
  if (alphaQuotientPtr != NULL) *alphaQuotientPtr = alphaQuotient;

  return std::min(1.,alphaQuotient);
}

// template<class V, class M>
// bool
// SequenceGenerator<V, M>::acceptAlpha(double alpha)
// {
//   bool result = false;
//
//   if      (alpha <= 0.                                ) result = false;
//   else if (alpha >= 1.                                ) result = true;
//   else if (alpha >= m_env.rngObject()->uniformSample()) result = true;
//   else                                                  result = false;
//
//   return result;
// }

template<class V, class M>
int
SequenceGenerator<V,M>::writeInfo(const BaseVectorSequence<V,M> & workingChain,
                                  std::ofstream & ofsvar) const
{
  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "\n"
                            << "\n-----------------------------------------------------"
                            << "\n Writing more information about the Markov chain " << workingChain.name() << " to output file ..."
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  int iRC = UQ_OK_RC;

  if (m_optionsObj->m_rawChainGenerateExtra) {
    // Write m_logTargets
    ofsvar << m_optionsObj->m_prefix << "logTargets_sub" << m_env.subIdString() << " = zeros(" << m_logTargets.size()
           << ","                    << 1
           << ");"
           << std::endl;
    ofsvar << m_optionsObj->m_prefix << "logTargets_sub" << m_env.subIdString() << " = [";
    for (unsigned int i = 0; i < m_logTargets.size(); ++i) {
      ofsvar << m_logTargets[i]
             << std::endl;
    }
    ofsvar << "];\n";

    // Write m_alphaQuotients
    ofsvar << m_optionsObj->m_prefix << "alphaQuotients_sub" << m_env.subIdString() << " = zeros(" << m_alphaQuotients.size()
           << ","                    << 1
           << ");"
           << std::endl;
    ofsvar << m_optionsObj->m_prefix << "alphaQuotients_sub" << m_env.subIdString() << " = [";
    for (unsigned int i = 0; i < m_alphaQuotients.size(); ++i) {
      ofsvar << m_alphaQuotients[i]
             << std::endl;
    }
    ofsvar << "];\n";
  }

  if (false) { // Don't see need for code below. Let it there though, compiling, in case it is needed in the future.
    // Write names of components
    ofsvar << m_optionsObj->m_prefix << "componentNames = {";
    m_vectorSpace.printComponentsNames(ofsvar,false);
    ofsvar << "};\n";
  }

  if ((m_env.subDisplayFile()                   ) &&
      (m_optionsObj->m_totallyMute == false)) {
    *m_env.subDisplayFile() << "\n-----------------------------------------------------"
                            << "\n Finished writing more information about the Markov chain " << workingChain.name()
                            << "\n-----------------------------------------------------"
                            << "\n"
                            << std::endl;
  }

  return iRC;
}

template <class V, class M>
const BaseTKGroup<V, M> &
SequenceGenerator<V, M>::transitionKernel() const
{
  return *m_tk;
}

template <class V, class M>
void
SequenceGenerator<V, M>::readFullChain(const std::string & inputFileName,
                                       const std::string & inputFileType,
                                       unsigned int chainSize,
                                       BaseVectorSequence<V, M> & workingChain)
{
  workingChain.unifiedReadContents(inputFileName, inputFileType, chainSize);
}

}  // End namespace QUESO

template class QUESO::SequenceGenerator<QUESO::GslVector, QUESO::GslMatrix>;
