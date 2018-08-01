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

#ifndef QUESO_ALGORITHM_FACTORY_H
#define QUESO_ALGORITHM_FACTORY_H

#include <queso/Factory.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/Algorithm.h>

namespace QUESO
{

/**
 * AlgorithmFactory class defintion.  Clients subclassing this for their own
 * algorithm (aka Metropolis-Hastings acceptance ratio) should implement
 * build_algorithm.
 */
template <typename V = GslVector, typename M = GslMatrix>
class AlgorithmFactory : public Factory<Algorithm<V, M> >
{
public:
  /**
   * Constructor. Takes the name to be mapped.
   */
  AlgorithmFactory(const std::string & name)
    : Factory<Algorithm<V, M> >(name)
  {}

  /**
   * Destructor. (Empty.)
   */
  virtual ~AlgorithmFactory() {}

  static void set_environment(const BaseEnvironment & env)
  {
    m_env = &env;
  }

  static void set_tk(const BaseTKGroup<V, M> & tk)
  {
    m_tk = &tk;
  }

protected:
  virtual typename SharedPtr<Algorithm<V, M> >::Type build_algorithm() = 0;

  static const BaseEnvironment * m_env;
  static const BaseTKGroup<V, M> * m_tk;

private:
  /**
   * Create a Base class.  Force this to be implemented
   * later.
   */
  virtual typename SharedPtr<Algorithm<V, M> >::Type create();
};

template <typename V, typename M>
inline
typename SharedPtr<Algorithm<V, M> >::Type
AlgorithmFactory<V, M>::create()
{
  queso_require_msg(m_env, "ERROR: must call set_environment() before building alg!");
  queso_require_msg(m_tk, "ERROR: must call set_tk() before building alg!");

  typename SharedPtr<Algorithm<V, M> >::Type new_alg = this->build_algorithm();

  queso_assert(new_alg);

  return new_alg;
}

/**
 * AlgorithmFactoryImp implementation of AlgorithmFactory.  Implements an
 * algorithm factory for the standard Metropolis-Hastings algorithm (aka
 * acceptance ratio).
 */
template <template <typename, typename> class DerivedAlgorithm, typename V = GslVector, typename M = GslMatrix>
class AlgorithmFactoryImp : public AlgorithmFactory<V, M>
{
public:
  AlgorithmFactoryImp(const std::string & name)
    : AlgorithmFactory<V, M>(name)
  {}

  virtual ~AlgorithmFactoryImp() {}

private:
  virtual typename SharedPtr<Algorithm<V, M> >::Type build_algorithm()
  {
    typename SharedPtr<Algorithm<V, M> >::Type new_alg;
    new_alg.reset(new DerivedAlgorithm<V, M>(*(this->m_env), *(this->m_tk)));
    return new_alg;
  }

};

} // namespace QUESO

#endif // QUESO_ALGORITHM_FACTORY_H
