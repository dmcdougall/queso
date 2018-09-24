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

#include "config_queso.h"

#ifdef QUESO_HAVE_CPPUNIT

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include <queso/EnvironmentOptions.h>
#include <queso/GslVector.h>
#include <queso/VectorSpace.h>
#include <queso/BoxSubset.h>
#include <queso/ConstantVectorFunction.h>
#include <queso/GslNumericVector.h>
#include <queso/GslSparseMatrix.h>

namespace QUESOTesting
{

class ConstantVectorFunctionTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(ConstantVectorFunctionTest);
  CPPUNIT_TEST(test_compute);
  CPPUNIT_TEST_SUITE_END();

  // yes, this is necessary
public:
  void setUp()
  {
    env.reset(new QUESO::FullEnvironment("","",NULL));

    space.reset(new QUESO::VectorSpace<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> >(*env, "", 1, NULL));

    min1.reset(new QUESO::GslNumericVector<libMesh::Number>(space->zeroVector()));
    (*min1)[0] = 1;
    max1.reset(new QUESO::GslNumericVector<libMesh::Number>(space->zeroVector()));
    (*max1)[0] = 3;

    domain.reset(new QUESO::BoxSubset<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> >("", *space, *min1, *max1));
    image.reset(new QUESO::BoxSubset<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> >("", *space, *min1, *max1));

    val.reset(new QUESO::GslNumericVector<libMesh::Number>(space->zeroVector()));
    (*val)[0] = 2.0;

    fn.reset(new QUESO::ConstantVectorFunction<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number>, QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> >("", *domain, *image, *val));
  }

  void test_compute()
  {
    QUESO::GslNumericVector<libMesh::Number> point(space->zeroVector());
    point[0] = 1.5;

    QUESO::GslNumericVector<libMesh::Number> result(space->zeroVector());
    fn->compute(point, NULL, result, NULL, NULL, NULL);

    CPPUNIT_ASSERT_EQUAL(result[0], (*val)[0]);
  }

private:
  typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type env;
  typename QUESO::ScopedPtr<QUESO::VectorSpace<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > >::Type space;
  typename QUESO::ScopedPtr<QUESO::GslNumericVector<libMesh::Number>>::Type min1;
  typename QUESO::ScopedPtr<QUESO::GslNumericVector<libMesh::Number>>::Type max1;
  typename QUESO::ScopedPtr<QUESO::BoxSubset<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > >::Type domain;
  typename QUESO::ScopedPtr<QUESO::BoxSubset<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > >::Type image;
  typename QUESO::ScopedPtr<QUESO::GslNumericVector<libMesh::Number>>::Type val;
  typename QUESO::ScopedPtr<QUESO::ConstantVectorFunction<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number>, QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > >::Type fn;
};

CPPUNIT_TEST_SUITE_REGISTRATION(ConstantVectorFunctionTest);

}  // end namespace QUESOTesting

#endif  // QUESO_HAVE_CPPUNIT
