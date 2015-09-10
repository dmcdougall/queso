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

#ifndef QUESO_SG_H
#define QUESO_SG_H

#include <iostream>
#include <vector>

namespace QUESO {

class GslVector;
class GslMatrix;
class BaseEnvironment;
class MhOptionsValues;
class MetropolisHastingsSGOptions;
template <class V, class M> class BaseTKGroup;
template <class V, class M> class BaseVectorRV;
template <class V, class M> class BaseVectorSequence;
template <class V, class M> class BaseJointPdf;
template <class V, class M> class VectorSpace;
template <class V, class M> class ScalarFunctionSynchronizer;
template <class V> class MarkovChainPositionData;
template <class T> class ScalarSequence;

//--------------------------------------------------
// MetropolisHastingsSG----------------------
//--------------------------------------------------

/*!\class MetropolisHastingsSG
 * \brief A templated class that represents a Metropolis-Hastings generator of samples.
 *
 * This class implements a Metropolis-Hastings generator of samples. 'SG' stands for 'Sequence Generator'.
 * Options reading is handled by class 'MetropolisHastingsOptions'. If options request data to be
 * written in the output file (MATLAB .m format only, for now), the user can check which MATLAB variables
 * are defined and set by running 'grep zeros <OUTPUT FILE NAME>' after the solution procedures ends.
 * The names of the variables are self explanatory. */

template <class V = GslVector, class M = GslMatrix>
class SequenceGenerator
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor.
  /*! This method reads the options from the options input file. It calls commonConstructor().
   * Requirements: 1) the image set of the vector random variable 'sourceRv' should belong to a
   * vector space of dimension equal to the size of the vector 'initialPosition' and 2) if
   * 'inputProposalCovMatrix' is not NULL, is should be square and its size should be equal to the
   * size of 'initialPosition'. If the requirements are satisfied, the constructor then reads input
   * options that begin with the string '\<prefix\>_mh_'. For instance, if 'prefix' is
   * 'pROblem_775_ip_', then the constructor will read all options that begin with 'pROblem_775_ip_mh_'.
    Options reading is handled by class 'MetropolisHastingsOptions'.*/
  SequenceGenerator(const char * prefix,
                    const MhOptionsValues * alternativeOptionsValues,
                    const BaseVectorRV<V, M> & sourceRv,
                    const V & initialPosition,
                    const M * inputProposalCovMatrix);

  //! Destructor
  virtual ~SequenceGenerator();
  //@}

  //! @name Statistical methods
  //@{
 //! Method to generate the chain.
 /*! Requirement: the vector space 'm_vectorSpace' should have dimension equal to the size of a
  * vector in 'workingChain'. If the requirement is satisfied, this operation sets the size and
  * the contents of 'workingChain' using the algorithm options set in the constructor. If not NULL,
  * 'workingLogLikelihoodValues' and 'workingLogTargetValues' are set accordingly. This operation
  * currently implements the DRAM algorithm (Heikki Haario, Marko Laine, Antonietta Mira and
  * Eero Saksman, "DRAM: Efficient Adaptive MCMC", Statistics and Computing (2006), 16:339-354),
  * as a translation of the core routine at the MCMC toolbox for MATLAB, available at
  *  \htmlonly helios.fmi.fi/~lainema/mcmc/‎ \endhtmlonly (accessed in July 3rd, 2013). Indeed,
  * the example available in examples/statisticalInverseProblem is closely related to the
  * 'normal example' in the toolbox. Over time, though:
  <list type=number>
  <item> the whole set of QUESO classes took shape, focusing not only on Markov Chains, but on
  statistical forward problems and model validation as well;
  <item> the interfaces to this Metropolis-Hastings class changed;
  <item> QUESO has parallel capabilities;
  <item> TK (transition kernel) class has been added in order to have both DRAM with Stochastic
  Newton capabilities.
  </list>
  */
  virtual void generateSequence(BaseVectorSequence<V, M> & workingChain,
                                ScalarSequence<double> * workingLogLikelihoodValues,
                                ScalarSequence<double> * workingLogTargetValues) = 0;

  //! Returns the underlying transition kernel for this sequence generator
  virtual const BaseTKGroup<V, M> & transitionKernel() const;
  //@}

  //! @name I/O methods
  //@{
  //! TODO: Prints the sequence.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
  friend std::ostream& operator<<(std::ostream& os,
      const SequenceGenerator<V, M>& obj)
  {
    obj.print(os);

    return os;
  }
  //@}

protected:
  //!  This method generates the chain.
  /*! Given the size of the chain, the values of the first position in the chain. It checks if the position
   * is out of target pdf support, otherwise it calculates the value of both its likelihood and prior and
   * adds the position to the chain. For the next positions, once they are generated, some tests are
   * performed (such as unicity and the value of alpha) and the steps for the first position are repeated,
   * including the optional Delayed Rejection and the adaptive Metropolis (adaptation of covariance matrix)
   * steps.*/
  virtual void generateFullChain(const V & valuesOf1stPosition,
                                 unsigned int chainSize,
                                 BaseVectorSequence<V, M> & workingChain,
                                 ScalarSequence<double> * workingLogLikelihoodValues,
                                 ScalarSequence<double> * workingLogTargetValues) = 0;

  //! This method reads the chain contents.
  void readFullChain(const std::string & inputFileName,
                     const std::string & inputFileType,
                     unsigned int chainSize,
                     BaseVectorSequence<V, M> & workingChain);

  //! Calculates acceptance ration.
  /*! It is called by alpha(const std::vector<MarkovChainPositionData<V>*>& inputPositions,
      const std::vector<unsigned int>& inputTKStageIds); */
  double alpha(const MarkovChainPositionData<V> & x,
               const MarkovChainPositionData<V> & y,
               unsigned int xStageId,
               unsigned int yStageId,
               double * alphaQuotientPtr = NULL);

  //! Decides whether or not to accept alpha.
  /*! If either alpha is negative or greater than one, its value will not be accepted.*/
  bool acceptAlpha(double alpha);

  //! Writes information about the Markov chain in a file.
  /*! It writes down the alpha quotients, number of positions out of
   * target support, the name of the components and the chain runtime.*/
  int writeInfo(const BaseVectorSequence<V, M> & workingChain,
                std::ofstream & ofsvar) const;

  const BaseEnvironment & m_env;
  const VectorSpace <V, M> & m_vectorSpace;
  const BaseJointPdf<V, M> & m_targetPdf;
  V m_initialPosition;
  M m_initialProposalCovMatrix;
  bool m_nullInputProposalCovMatrix;
  const ScalarFunctionSynchronizer<V, M> * m_targetPdfSynchronizer;
  BaseTKGroup<V, M> * m_tk;
  std::vector<unsigned int> m_idsOfUniquePositions;
  std::vector<double> m_logTargets;
  std::vector<double> m_alphaQuotients;
  double m_lastChainSize;
  V * m_lastMean;
  M * m_lastAdaptedCovMatrix;
  const MhOptionsValues * m_optionsObj;
  MetropolisHastingsSGOptions * m_oldOptions;
  bool m_computeInitialPriorAndLikelihoodValues;
  double m_initialLogPriorValue;
  double m_initialLogLikelihoodValue;
  bool m_userDidNotProvideOptions;

private:
  //! Reads the options values from the options input file.
  /*!  This method \b actually reads the options input file, such as the value for the Delayed Rejection
   * scales, the presence of Hessian covariance matrices and reads the user-provided initial proposal
   * covariance matrix (alternative to local Hessians).*/
  void commonConstructor();

  //! Adaptive Metropolis method that deals with adapting the proposal covariance matrix
  virtual void adapt(unsigned int positionId,
                     BaseVectorSequence<V, M> & workingChain) = 0;

  //! This method updates the adapted covariance matrix
  /*! This function is called is the option to used adaptive Metropolis was chosen by the user
   * (via options input file). It performs an adaptation of covariance matrix. */
  virtual void updateAdaptedCovMatrix(const BaseVectorSequence<V, M> & subChain,
                                      unsigned int idOfFirstPositionInSubChain,
                                      double & lastChainSize,
                                      V & lastMean,
                                      M & lastAdaptedCovMatrix) = 0;
};

}  // End namespace QUESO

#endif  // QUESO_SG_H
