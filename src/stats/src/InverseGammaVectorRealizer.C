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

#include <limits>
#include <queso/InverseGammaVectorRealizer.h>
#include <queso/GslVector.h>
#include <queso/GslNumericVector.h>
#include <queso/GslMatrix.h>
#include <queso/GslSparseMatrix.h>
#include <libmesh/eigen_sparse_vector.h>
#include <libmesh/eigen_sparse_matrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V, class M>
InverseGammaVectorRealizer<V,M>::InverseGammaVectorRealizer(
  const char*                  prefix,
  const VectorSet<V,M>& unifiedImageSet,
  const V&                     alpha,
  const V&                     beta)
  :
  BaseVectorRealizer<V,M>(((std::string)(prefix)+"gen").c_str(),unifiedImageSet,std::numeric_limits<unsigned int>::max()),
  m_alpha(alpha),
  m_beta (beta)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering InverseGammaVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving InverseGammaVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
InverseGammaVectorRealizer<V,M>::~InverseGammaVectorRealizer()
{
}
// Realization-related methods----------------------
template<class V, class M>
void
InverseGammaVectorRealizer<V,M>::realization(V& nextValues) const
{
  nextValues.cwSetInverseGamma(m_alpha,m_beta);
  return;
}

template class InverseGammaVectorRealizer<GslVector, GslMatrix>;
template class InverseGammaVectorRealizer<GslNumericVector<libMesh::Number>, GslSparseMatrix<libMesh::Number> >;
template class InverseGammaVectorRealizer<libMesh::EigenSparseVector<libMesh::Number>, libMesh::EigenSparseMatrix<libMesh::Number> >;

}  // End namespace QUESO
