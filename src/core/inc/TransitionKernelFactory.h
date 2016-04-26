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

#ifndef QUESO_TK_FACTORY_H
#define QUESO_TK_FACTORY_H

#include <queso/FactoryWithVectorSpace.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/TKGroup.h>
#include <queso/MetropolisHastingsSGOptions.h>

namespace QUESO
{

/**
 * TransitionKernelFactory class defintion.
 */
class TransitionKernelFactory : public FactoryWithVectorSpace<BaseTKGroup<GslVector, GslMatrix> >
{
public:
  /**
   * Constructor. Takes the name to be mapped.
   */
  TransitionKernelFactory(const std::string & name)
    : FactoryWithVectorSpace<BaseTKGroup<GslVector, GslMatrix> >(name)
  {}

  /**
   * Destructor. (Empty.)
   */
  virtual ~TransitionKernelFactory() {}

protected:
  static void set_dr_scales(const std::vector<double> & scales)
  {
    m_dr_scales = &scales;
  }

  static void set_pdf_synchronizer(const ScalarFunctionSynchronizer<GslVector, GslMatrix> & synchronizer)
  {
    m_pdf_synchronizer = &synchronizer;
  }

  static void set_initial_cov_matrix(const GslMatrix & cov_matrix)
  {
    m_initial_cov_matrix = &cov_matrix;
  }

  static void set_options(const MhOptionsValues & options)
  {
    m_options = &options;
  }

  virtual SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type build_tk(
      const MhOptionsValues & options,
      const VectorSpace<GslVector, GslMatrix> & v,
      const std::vector<double> & dr_scales,
      const ScalarFunctionSynchronizer<GslVector, GslMatrix> & pdf_synchronizer,
      const GslMatrix & initial_cov_matrix);

  static const std::vector<double> * m_dr_scales;
  static const ScalarFunctionSynchronizer<GslVector, GslMatrix> * m_pdf_synchronizer;
  static const GslMatrix * m_initial_cov_matrix;
  static const MhOptionsValues * m_options;

private:
  /**
   * Create a Base class.  Force this to be implemented
   * later.
   */
  virtual typename SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type create();
};

inline
typename SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type TransitionKernelFactory::create()
{
  queso_require_msg(m_vectorSpace, "ERROR: must call set_vectorspace() before building tk!");
  queso_require_msg(m_dr_scales, "ERROR: must call set_dr_scales() before building tk!");
  queso_require_msg(m_pdf_synchronizer, "ERROR: must call set_pdf_synchronizer() before building tk!");
  queso_require_msg(m_initial_cov_matrix, "ERROR: must call set_initial_cov_matrix() before building tk!");
  queso_require_msg(m_options, "ERROR: must call set_options() before building tk!");

  SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type new_tk = this->build_tk(
      *m_options,
      *m_vectorSpace,
      *m_dr_scales,
      *m_pdf_synchronizer,
      *m_initial_cov_matrix);

  queso_assert(new_tk);

  return new_tk;
}

} // namespace QUESO

#endif // QUESO_TK_FACTORY_H
