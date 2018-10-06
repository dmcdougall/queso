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
#ifndef EX_QOI_H
#define EX_QOI_H

#include <queso/GslMatrix.h>
#include <queso/DistArray.h>
#include <queso/GslNumericVector.h>
#include <queso/GslSparseMatrix.h>

struct
qoiRoutine_DataType
{
  double coef1;
  double coef2;
};

void
qoiRoutine(
  const QUESO::GslNumericVector<libMesh::Number>&                    paramValues,
  const QUESO::GslNumericVector<libMesh::Number>*                    paramDirection,
  const void*                                functionDataPtr,
        QUESO::GslNumericVector<libMesh::Number>&                    qoiValues,
        QUESO::DistArray<QUESO::GslNumericVector<libMesh::Number>*>* gradVectors,
        QUESO::DistArray<QUESO::GslSparseMatrix<libMesh::Number>*>* hessianMatrices,
        QUESO::DistArray<QUESO::GslNumericVector<libMesh::Number>*>* hessianEffects);

#endif
