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
class MpiComm;
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
// MHRawChainInfoStruct --------------------------
//--------------------------------------------------
 /*! \file MetropolisHastingsSG.h
 * \brief A templated class that represents a Metropolis-Hastings generator of samples and a struct which stores its info.
 *
 * \struct  MHRawChainInfoStruct
 * \brief A struct that represents a Metropolis-Hastings sample.
 *
 * Some of the information about the  Metropolis-Hastings sample generator includes the allowed number
 * of delayed rejections, number of rejections, number of positions in or out of target support, and
 * so on. This struct is responsible for the storage of such info. */

struct MHRawChainInfoStruct
{
 //! @name Constructor/Destructor methods
 //@{
 //! Constructor.
  MHRawChainInfoStruct();

  //! Copy constructor.
  MHRawChainInfoStruct(const MHRawChainInfoStruct& rhs);

  //! Destructor
  ~MHRawChainInfoStruct();
  //@}

  //! @name Set methods
  //@{
  //! Assignment operator.
  MHRawChainInfoStruct& operator= (const MHRawChainInfoStruct& rhs);

  //! Addition assignment operator.
  MHRawChainInfoStruct& operator+=(const MHRawChainInfoStruct& rhs);
  //@}

   //! @name Misc methods
  //@{
  //! Copies Metropolis-Hastings chain info from \c src to \c this.
  void copy  (const MHRawChainInfoStruct& src);

  //! Resets Metropolis-Hastings chain info.
  void reset ();

  //! Calculates the MPI sum of \c this.
  void mpiSum(const MpiComm& comm, MHRawChainInfoStruct& sumInfo);
  //@}

  double       runTime;
  double       candidateRunTime;
  double       targetRunTime;
  double       mhAlphaRunTime;
  double       drAlphaRunTime;
  double       drRunTime;
  double       amRunTime;

  unsigned int numTargetCalls;
  unsigned int numDRs;
  unsigned int numOutOfTargetSupport;
  unsigned int numOutOfTargetSupportInDR;
  unsigned int numRejections;

};

/*!
 * \class SequenceGenerator
 * \brief A templated class that represents an abstract base sequence generator
 * of samples.
 *
 * This class provides an interface for sequence generators (samplers).
 * Options reading is handled by class 'MhOptionsValues'.
 */
template <class V = GslVector, class M = GslMatrix>
class SequenceGenerator
{
public:
  //!@name Constructor/Destructor methods
  //@{
  //! Constructor.
  /*!
   * This method reads the options from the options input file.  Requirements:
   * -# the image set of the vector random variable 'sourceRv' should belong to
   *    a vector space of dimension equal to the size of the vector
   *    'initialPosition', and;
   * -# if 'inputProposalCovMatrix' is not NULL, is should be square and its
   *    size should be equal to the size of 'initialPosition'.
   * .
   *
   * If the requirements are satisfied, the constructor then reads input
   * options that begin with the string '\<prefix\>_mh_'. For instance, if
   * 'prefix' is 'pROblem_775_ip_', then the constructor will read all options
   * that begin with 'pROblem_775_ip_mh_'.  Options reading is handled by class
   * 'MhOptionsValues'.
   */
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
  /*!
   * Requirement: the vector space 'm_vectorSpace' should have dimension equal to the size of a
   * vector in 'workingChain'. If the requirement is satisfied, this operation sets the size and
   * the contents of 'workingChain' using the algorithm options set in the constructor. If not NULL,
   * 'workingLogLikelihoodValues' and 'workingLogTargetValues' are set accordingly. This operation
   * currently implements the DRAM algorithm (Heikki Haario, Marko Laine, Antonietta Mira and
   * Eero Saksman, "DRAM: Efficient Adaptive MCMC", Statistics and Computing (2006), 16:339-354),
   * as a translation of the core routine at the MCMC toolbox for MATLAB, available at
   *  \htmlonly helios.fmi.fi/~lainema/mcmc/‎ \endhtmlonly (accessed in July 3rd, 2013). Indeed,
   * the example available in examples/statisticalInverseProblem is closely related to the
   * 'normal example' in the toolbox. Over time, though:
   * <list type=number>
   * <item> the whole set of QUESO classes took shape, focusing not only on Markov Chains, but on
   * statistical forward problems and model validation as well;
   * <item> the interfaces to this Metropolis-Hastings class changed;
   * <item> QUESO has parallel capabilities;
   * <item> TK (transition kernel) class has been added in order to have both DRAM with Stochastic
   * Newton capabilities.
   * </list>
   *
   * This method is responsible for setting the size of the \c workingChain,
   * \c workingLogTargetValues, and \c workingLogTargetValues
   */
  virtual void generateSequence(BaseVectorSequence<V, M> & workingChain,
                                ScalarSequence<double> * workingLogLikelihoodValues,
                                ScalarSequence<double> * workingLogTargetValues);

  void generateFullChain(const V & valuesOf1stPosition,
                         unsigned int chainSize,
                         BaseVectorSequence<V, M> & workingChain,
                         ScalarSequence<double> * workingLogLikelihoodValues,
                         ScalarSequence<double> * workingLogTargetValues);

  //! Implement this method to compute the proposed state of the chain
  /*!
   * \c positionId is the index of the current state, and \c
   * workingChain is the whole chain.  Therefore the index of the proposed
   * state is \c positionId + 1.
   *
   * The reason the whole working chain is passed is because one may wish to
   * use chain history to compute the proposed state (cf. adaptive Metropolis)
   * and we do not want to disallow this.
   *
   * The implementation may should set \c proposedState so QUESO may then
   * decide whether or not to accept this sample by the Metropolis-Hastings
   * algorithm.
   */
  virtual void propose(unsigned int positionId,
                       const BaseVectorSequence<V, M> & workingChain,
                       V & proposedState) = 0;

  //! Do delayed rejection
  virtual bool delayedRejection(unsigned int positionId,
      MarkovChainPositionData<V> & currentPositionData,
      MarkovChainPositionData<V> & currentCandidateData) = 0;

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

  //! Calculates acceptance ration.
  /*! The acceptance ratio is used to decide whether to accept or reject a candidate. */
  virtual double alpha(
      const std::vector<MarkovChainPositionData<V>*>& inputPositions,
      const std::vector<unsigned int>& inputTKStageIds);

  //! Decides whether or not to accept alpha.
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
};

}  // End namespace QUESO

#endif  // QUESO_SG_H
