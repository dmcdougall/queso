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

// This class
#include <queso/InterpolationSurrogateIOBase.h>

// QUESO
#include <queso/VectorSpace.h>
#include <queso/GslVector.h>
#include <queso/GslNumericVector.h>
#include <queso/GslMatrix.h>
#include <queso/GslSparseMatrix.h>
#include <libmesh/eigen_sparse_vector.h>
#include <libmesh/eigen_sparse_matrix.h>

namespace QUESO
{
  template<class V, class M>
  InterpolationSurrogateIOBase<V,M>::InterpolationSurrogateIOBase()
  {}

  // Instantiate
  template class InterpolationSurrogateIOBase<GslVector,GslMatrix>;
  template class InterpolationSurrogateIOBase<GslNumericVector<libMesh::Number>,GslSparseMatrix<libMesh::Number> >;
  template class InterpolationSurrogateIOBase<libMesh::EigenSparseVector<libMesh::Number>,libMesh::EigenSparseMatrix<libMesh::Number> >;

} // end namespace QUESO
