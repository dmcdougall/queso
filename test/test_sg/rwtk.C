#include "rwtk.h"

#include <queso/VectorSpace.h>
#include <queso/VectorRV.h>

template <class V, class M>
RWTK<V, M>::RWTK(const char * prefix,
                 const QUESO::VectorSpace<V, M> & vectorSpace,
                 const std::vector<double> & scales)
: QUESO::BaseTKGroup<V, M>(prefix, vectorSpace, scales)
{
  // Do nothing?
}

template <class V, class M>
RWTK<V, M>::~RWTK()
{
  // Do nothing?
}

template <class V, class M>
bool
RWTK<V, M>::symmetric() const
{
  return false;
}

template <class V, class M>
const QUESO::BaseVectorRV<V, M> &
RWTK<V, M>::rv(unsigned int stageId) const
{
  return *(this->m_rv);
}

template <class V, class M>
const QUESO::BaseVectorRV<V, M> &
RWTK<V, M>::rv(const std::vector<unsigned int> & stageIds)
{
  return *(this->m_rv);
}

template class RWTK<QUESO::GslVector, QUESO::GslMatrix>;
