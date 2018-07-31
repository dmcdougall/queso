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

#ifndef QUESO_TK_FACTORY_H
#define QUESO_TK_FACTORY_H

#include <queso/Factory.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/TKGroup.h>
#include <queso/MetropolisHastingsSGOptions.h>

namespace QUESO
{

/**
 * TransitionKernelFactory class defintion.
 */
template <typename V = GslVector, typename M = GslMatrix>
class TransitionKernelFactory : public Factory<BaseTKGroup<V, M> >
{
public:
  /**
   * Constructor. Takes the name to be mapped.
   */
  TransitionKernelFactory(const std::string & name)
    : Factory<BaseTKGroup<V, M> >(name)
  {}

  /**
   * Destructor. (Empty.)
   */
  virtual ~TransitionKernelFactory() {}

  /**
   * Static method to set the vector space the transition kernel is defined on
   */
  static void set_vectorspace(const VectorSpace<V, M> & v)
  {
    m_vectorSpace = &v;
  }

  /**
   * Static method to set the vector of scale factor to scale a proposal
   * covariance matrix by for the purpose of delayed rejection
   */
  static void set_dr_scales(const std::vector<double> & scales)
  {
    m_dr_scales = &scales;
  }

  /**
   * Static method to set the pdf synchronizer.  Used by Stochastic Newton?
   */
  static void set_pdf_synchronizer(const ScalarFunctionSynchronizer<V, M> & synchronizer)
  {
    m_pdf_synchronizer = &synchronizer;
  }

  /**
   * Static method to set the initial proposal covariance matrix
   */
  static void set_initial_cov_matrix(M & cov_matrix)
  {
    m_initial_cov_matrix = &cov_matrix;
  }

  /**
   * Static method to set the options object.  Useful for passing input file
   * options to transition kernels.
   */
  static void set_options(const MhOptionsValues & options)
  {
    m_options = &options;
  }

  /**
   * Static method to set the pdf we wish to draw samples from
   */
  static void set_target_pdf(const BaseJointPdf<V, M> & target_pdf)
  {
    m_target_pdf = &target_pdf;
  }

protected:
  virtual typename SharedPtr<BaseTKGroup<V, M> >::Type build_tk() = 0;

  static const VectorSpace<V, M> * m_vectorSpace;
  static const std::vector<double> * m_dr_scales;
  static const ScalarFunctionSynchronizer<V, M> * m_pdf_synchronizer;
  static M * m_initial_cov_matrix;
  static const MhOptionsValues * m_options;
  static const BaseJointPdf<V, M> * m_target_pdf;

private:
  /**
   * Create a Base class.  Force this to be implemented
   * later.
   */
  virtual typename SharedPtr<BaseTKGroup<V, M> >::Type create();
};

template <typename V, typename M>
inline
typename SharedPtr<BaseTKGroup<V, M> >::Type TransitionKernelFactory<V, M>::create()
{
  queso_require_msg(m_vectorSpace, "ERROR: must call set_vectorspace() before building tk!");
  queso_require_msg(m_dr_scales, "ERROR: must call set_dr_scales() before building tk!");
  queso_require_msg(m_pdf_synchronizer, "ERROR: must call set_pdf_synchronizer() before building tk!");
  queso_require_msg(m_initial_cov_matrix, "ERROR: must call set_initial_cov_matrix() before building tk!");
  queso_require_msg(m_options, "ERROR: must call set_options() before building tk!");
  queso_require_msg(m_target_pdf, "ERROR: must call set_target_pdf() before building tk!");

  typename SharedPtr<BaseTKGroup<V, M> >::Type new_tk = this->build_tk();

  queso_assert(new_tk);

  return new_tk;
}

} // namespace QUESO

#endif // QUESO_TK_FACTORY_H
