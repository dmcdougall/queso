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

#include <queso/Environment.h>
#include <queso/ScopedPtr.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/SequenceOfVectors.h>
#include <queso/BoxSubset.h>
#include <queso/GslNumericVector.h>
#include <queso/GslSparseMatrix.h>

#include <cmath>
#include <vector>
#include <string>

namespace QUESOTesting
{

class SequenceOfVectorsTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(SequenceOfVectorsTest);
  CPPUNIT_TEST(test_uniformly_sampled_cdf);
  CPPUNIT_TEST(test_std_var);
  CPPUNIT_TEST(test_hists);
  CPPUNIT_TEST(test_iqr);
  CPPUNIT_TEST(test_auto_covariance);
  CPPUNIT_TEST(test_read);
  CPPUNIT_TEST(test_scale_kde);
  CPPUNIT_TEST(test_gaussian_kde);
  CPPUNIT_TEST_SUITE_END();

  // yes, this is necessary
public:
  void setUp()
  {
    env.reset(new QUESO::FullEnvironment("","",NULL));
    space.reset(new QUESO::VectorSpace<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> >(*env, "", 2, NULL));
    sequence.reset(new QUESO::SequenceOfVectors<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> >(*space, 13, ""));

    // Fill up the sequence with some not-so-random stuff
    QUESO::GslNumericVector<libMesh::Number> v(space->zeroVector());
    for (unsigned int i = 0; i < sequence->subSequenceSize(); i++) {
      v[0] = i;
      v[1] = i + 1;
      sequence->setPositionValues(i, v);
    }
  }

  void test_uniformly_sampled_cdf()
  {
    QUESO::GslNumericVector<libMesh::Number> num_points(space->zeroVector());
    num_points.cwSet(14.0);

    QUESO::VectorSpace<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > rowSpace(*env, "", sequence->subSequenceSize(), NULL);
    QUESO::ArrayOfOneDGrids<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > cdfGrids("", rowSpace);
    QUESO::ArrayOfOneDTables<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > cdfValues("", rowSpace);
    sequence->subUniformlySampledCdf(num_points, cdfGrids, cdfValues);

    for (unsigned int i = 0; i < sequence->subSequenceSize(); i++) {
      unsigned int cdf_val1 = cdfValues.oneDTable(0)[i] * 13;
      unsigned int cdf_val2 = cdfValues.oneDTable(1)[i] * 13;
      CPPUNIT_ASSERT_EQUAL(i, cdf_val1);
      CPPUNIT_ASSERT_EQUAL(i, cdf_val2);
    }

    sequence->unifiedUniformlySampledCdf(num_points, cdfGrids, cdfValues);

    for (unsigned int i = 0; i < sequence->subSequenceSize(); i++) {
      unsigned int cdf_val1 = cdfValues.oneDTable(0)[i] * 13;
      unsigned int cdf_val2 = cdfValues.oneDTable(1)[i] * 13;
      CPPUNIT_ASSERT_EQUAL(i, cdf_val1);
      CPPUNIT_ASSERT_EQUAL(i, cdf_val2);
    }
  }

  void test_std_var()
  {
    QUESO::GslNumericVector<libMesh::Number> mean(space->zeroVector());
    mean[0] = 6.0;
    mean[1] = 7.0;
    QUESO::GslNumericVector<libMesh::Number> std(space->zeroVector());
    QUESO::GslNumericVector<libMesh::Number> var(space->zeroVector());

    double actualSampleStd = std::sqrt(182.0 / 12.0);
    sequence->subSampleStd(0, sequence->subSequenceSize(), mean, std);

    CPPUNIT_ASSERT_EQUAL(actualSampleStd, std[0]);
    CPPUNIT_ASSERT_EQUAL(actualSampleStd, std[1]);

    sequence->unifiedSampleStd(0, sequence->subSequenceSize(), mean, std);

    CPPUNIT_ASSERT_EQUAL(actualSampleStd, std[0]);
    CPPUNIT_ASSERT_EQUAL(actualSampleStd, std[1]);

    double actualPopulationVar = 182.0 / 13.0;
    sequence->subPopulationVariance(0, sequence->subSequenceSize(), mean, var);

    CPPUNIT_ASSERT_EQUAL(actualPopulationVar, var[0]);
    CPPUNIT_ASSERT_EQUAL(actualPopulationVar, var[1]);

    sequence->unifiedPopulationVariance(0, sequence->subSequenceSize(), mean, var);

    CPPUNIT_ASSERT_EQUAL(actualPopulationVar, var[0]);
    CPPUNIT_ASSERT_EQUAL(actualPopulationVar, var[1]);
  }

  void test_hists()
  {
    QUESO::GslNumericVector<libMesh::Number> min(space->zeroVector());
    QUESO::GslNumericVector<libMesh::Number> max(space->zeroVector());
    std::vector<QUESO::GslNumericVector<libMesh::Number> *> centers(14, (QUESO::GslNumericVector<libMesh::Number> *)NULL);
    std::vector<QUESO::GslNumericVector<libMesh::Number> *> bins1(14, (QUESO::GslNumericVector<libMesh::Number> *)NULL);
    std::vector<QUESO::GslNumericVector<libMesh::Number> *> bins2(14, (QUESO::GslNumericVector<libMesh::Number> *)NULL);

    sequence->subMinMaxExtra(0, sequence->subSequenceSize(), min, max);
    sequence->subHistogram(0, min, max, centers, bins1);
    sequence->unifiedHistogram(0, min, max, centers, bins2);

    for (unsigned int i = 0; i < 14; i++) {
      CPPUNIT_ASSERT_EQUAL((*bins1[i])[0], (*bins2[i])[0]);
      CPPUNIT_ASSERT_EQUAL((*bins1[i])[1], (*bins2[i])[1]);
    }

    // The first element appears to always be zero.
    for (unsigned int i = 1; i < 14; i++) {
      CPPUNIT_ASSERT_EQUAL(1.0, (*bins1[i])[0]);
      CPPUNIT_ASSERT_EQUAL(1.0, (*bins1[i])[1]);
      CPPUNIT_ASSERT_EQUAL(1.0, (*bins2[i])[0]);
      CPPUNIT_ASSERT_EQUAL(1.0, (*bins2[i])[1]);
    }
  }

  void test_iqr()
  {
    QUESO::GslNumericVector<libMesh::Number> iqr(space->zeroVector());
    sequence->subInterQuantileRange(0, iqr);

    double actualIqr = 7.0;

    CPPUNIT_ASSERT_EQUAL(actualIqr, iqr[0]);
    CPPUNIT_ASSERT_EQUAL(actualIqr, iqr[1]);

    sequence->unifiedInterQuantileRange(0, iqr);

    CPPUNIT_ASSERT_EQUAL(actualIqr, iqr[0]);
    CPPUNIT_ASSERT_EQUAL(actualIqr, iqr[1]);
  }

  void test_auto_covariance()
  {
    QUESO::GslNumericVector<libMesh::Number> actualMean(space->zeroVector());
    actualMean[0] = 6.0;
    actualMean[1] = 7.0;
    double actualPopVar = 182.0 / 13.0;
    unsigned int numPos = sequence->subSequenceSize();

    // Compute auto covariance with zero lag.  Should be population var?
    QUESO::GslNumericVector<libMesh::Number> autoCov(space->zeroVector());
    sequence->autoCovariance(0, numPos, actualMean, 0, autoCov);
    CPPUNIT_ASSERT_EQUAL(actualPopVar, autoCov[0]);
    CPPUNIT_ASSERT_EQUAL(actualPopVar, autoCov[1]);

    // This is normalised by autoCov, so result should be 1.0 with lag 0
    QUESO::GslNumericVector<libMesh::Number> autoCorrViaDef(space->zeroVector());
    sequence->autoCorrViaDef(0, numPos, 0, autoCorrViaDef);
    CPPUNIT_ASSERT_EQUAL(1.0, autoCorrViaDef[0]);
    CPPUNIT_ASSERT_EQUAL(1.0, autoCorrViaDef[1]);

    std::vector<unsigned int> lags(2, 0);
    std::vector<QUESO::GslNumericVector<libMesh::Number> *> autoCorrsViaFft(2, (QUESO::GslNumericVector<libMesh::Number> *)NULL);
    sequence->autoCorrViaFft(0, numPos, lags, autoCorrsViaFft);
    CPPUNIT_ASSERT_EQUAL(1.0, (*autoCorrsViaFft[0])[0]);
    CPPUNIT_ASSERT_EQUAL(1.0, (*autoCorrsViaFft[0])[1]);

    QUESO::GslNumericVector<libMesh::Number> autoCorrsViaFftSum(space->zeroVector());
    sequence->autoCorrViaFft(0, numPos, 1, autoCorrsViaFftSum);
    CPPUNIT_ASSERT_EQUAL(1.0, autoCorrsViaFftSum[0]);
    CPPUNIT_ASSERT_EQUAL(1.0, autoCorrsViaFftSum[1]);
  }

  void test_read()
  {
    QUESO::SequenceOfVectors<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > read_sequence(*space, 5, "");

    // Filename logic
    std::string sequenceFileName = "unit/read_vector_sequence";
    const char * test_srcdir = std::getenv("srcdir");
    if (test_srcdir) {
      sequenceFileName = test_srcdir + ('/' + sequenceFileName);
    }

    read_sequence.unifiedReadContents(sequenceFileName, "m", 5);

    QUESO::GslNumericVector<libMesh::Number> value(space->zeroVector());
    for (unsigned int i = 0; i < 5; i++) {
      double val1 = i + 1;
      double val2 = i + 2;
      read_sequence.getPositionValues(i, value);
      CPPUNIT_ASSERT_EQUAL(val1, value[0]);
      CPPUNIT_ASSERT_EQUAL(val2, value[1]);
    }
  }

  void test_scale_kde()
  {
    double TOL = 1e-12;
    QUESO::GslNumericVector<libMesh::Number> actualIrq(space->zeroVector());
    actualIrq.cwSet(7.0);

    // Regression solution
    double actualKDEScale = 2.471509395452762269940194528317;

    QUESO::GslNumericVector<libMesh::Number> KDEScale(space->zeroVector());
    sequence->subScalesForKde(0, actualIrq, 1, KDEScale);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualKDEScale, KDEScale[0], TOL);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualKDEScale, KDEScale[1], TOL);

    sequence->unifiedScalesForKde(0, actualIrq, 1, KDEScale);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualKDEScale, KDEScale[0], TOL);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualKDEScale, KDEScale[1], TOL);
  }

  void test_gaussian_kde()
  {
    double TOL = 1e-12;

    QUESO::SequenceOfVectors<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > point(*space, 1, "");
    QUESO::GslNumericVector<libMesh::Number> mp1(space->zeroVector());
    mp1[0] = -0.1;
    mp1[1] = -0.1;
    point.setPositionValues(0, mp1);

    std::vector<QUESO::GslNumericVector<libMesh::Number> *> positions(1, &mp1);
    std::vector<QUESO::GslNumericVector<libMesh::Number> *> density(1, (QUESO::GslNumericVector<libMesh::Number> *)NULL);
    QUESO::GslNumericVector<libMesh::Number> scale(space->zeroVector());
    scale.cwSet(1.0);
    point.subGaussian1dKde(0, scale, positions, density);

    double actualKDE = 1.0 / std::sqrt(2.0 * M_PI);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualKDE, (*density[0])[0], TOL);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualKDE, (*density[0])[1], TOL);

    point.unifiedGaussian1dKde(0, scale, positions, density);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualKDE, (*density[0])[0], TOL);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(actualKDE, (*density[0])[1], TOL);
  }

private:
  typename QUESO::ScopedPtr<QUESO::BaseEnvironment>::Type env;
  typename QUESO::ScopedPtr<QUESO::VectorSpace<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > >::Type space;
  typename QUESO::ScopedPtr<QUESO::SequenceOfVectors<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > >::Type sequence;
};

CPPUNIT_TEST_SUITE_REGISTRATION(SequenceOfVectorsTest);

}  // end namespace QUESOTesting

#endif  // QUESO_HAVE_CPPUNIT
