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

#include <queso/AcceptanceProbability.h>

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/TKGroup.h>
#include <queso/ScalarFunction.h>
#include <queso/JointPdf.h>

namespace QUESO {

template <class V, class M>
AcceptanceProbability<V, M>::AcceptanceProbability(
    const BaseTKGroup<V, M> & transition_kernel,
    const BaseJointPdf<V, M> & prior_pdf,
    const BaseScalarFunction<V, M> & likelihood_fn)
:
  m_transition_kernel(transition_kernel),
  m_prior_pdf(prior_pdf),
  m_likelihood_fn(likelihood_fn)
{
}

template <class V, class M>
AcceptanceProbability<V, M>::~AcceptanceProbability()
{
}

template <class V, class M>
const BaseTKGroup<V, M> &
AcceptanceProbability<V, M>::transition_kernel() const
{
  return m_transition_kernel;
}

template <class V, class M>
const BaseJointPdf<V, M> &
AcceptanceProbability<V, M>::prior_pdf() const
{
  return m_prior_pdf;
}

template <class V, class M>
const BaseScalarFunction<V, M> &
AcceptanceProbability<V, M>::likelihood_function() const
{
  return m_likelihood_fn;
}

template class AcceptanceProbability<GslVector, GslMatrix>;
}  // End namespace QUESO
