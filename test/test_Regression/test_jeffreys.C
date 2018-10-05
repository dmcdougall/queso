//----------------------------------------------------------------
//----------------------------------------------------------------
//this is a test of the jeffreys sampling code
//----------------------------------------------------------------
//----------------------------------------------------------------

#include <queso/Environment.h>
#include <queso/GslMatrix.h>
#include <queso/StatisticalForwardProblem.h>
#include <queso/GenericVectorFunction.h>
#include <queso/JeffreysVectorRV.h>
#include <queso/GslMatrix.h>
#include <queso/DistArray.h>
#include <queso/GslNumericVector.h>
#include <queso/GslSparseMatrix.h>

#include <cstdlib>

struct qoiRoutine_DataType
{
  double coef1;
};

void qoiRoutine(
    const QUESO::GslNumericVector<libMesh::Number>&                    paramValues,
    const QUESO::GslNumericVector<libMesh::Number>*                    paramDirection,
    const void*                                functionDataPtr,
          QUESO::GslNumericVector<libMesh::Number>&                    qoiValues,
          QUESO::DistArray<QUESO::GslNumericVector<libMesh::Number>*>* gradVectors,
          QUESO::DistArray<QUESO::GslSparseMatrix<libMesh::Number>*>* hessianMatrices,
          QUESO::DistArray<QUESO::GslNumericVector<libMesh::Number>*>* hessianEffects)
{
  //logic to avoid warnings from intel compiler
  const QUESO::GslNumericVector<libMesh::Number>* aux1 = paramDirection;
  if (aux1) {};
  QUESO::DistArray<QUESO::GslNumericVector<libMesh::Number>*>* aux2 = gradVectors;
  if (aux2) {};
  aux2 = hessianEffects;
  QUESO::DistArray<QUESO::GslSparseMatrix<libMesh::Number>*>* aux3 = hessianMatrices;
  if (aux3) {};

  //checking size of paramValues and qoiValues
  //both should be size 1
  UQ_FATAL_TEST_MACRO(paramValues.sizeGlobal() != 1,
      QUESO::UQ_UNAVAILABLE_RANK,
      "qoiRoutine()",
      "paramValues vector does not have size 1");

  UQ_FATAL_TEST_MACRO(qoiValues.sizeGlobal() != 1,
      QUESO::UQ_UNAVAILABLE_RANK,
      "qoiRoutine()",
      "qoiValues vector does not have size 1");

 //compute qoi
 const QUESO::BaseEnvironment& env = paramValues.env();
 if(env.subRank() == 0) {
   double coef1 = ((qoiRoutine_DataType *) functionDataPtr)->coef1;
   qoiValues[0] = coef1*paramValues[0];
 }
 else {
   qoiValues[0] = 0.;
 }

 return;
}

void compute(const QUESO::FullEnvironment& env) {

  //step 1: instatiate parameter space
  QUESO::VectorSpace<QUESO::GslNumericVector<libMesh::Number>,QUESO::GslSparseMatrix<libMesh::Number>>
    paramSpace(env, "param_", 1, NULL);

  //step 2: instantiate the parameter domain
  QUESO::GslNumericVector<libMesh::Number> paramMins(paramSpace.zeroVector());
  paramMins.cwSet(0.001);
  QUESO::GslNumericVector<libMesh::Number> paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet(100.); //TODO: this is not working with gsl sampling right now in FP
  QUESO::BoxSubset<QUESO::GslNumericVector<libMesh::Number>,QUESO::GslSparseMatrix<libMesh::Number>>
    paramDomain("param_", paramSpace, paramMins, paramMaxs);

  //step 3: instantiate the qoi space
  QUESO::VectorSpace<QUESO::GslNumericVector<libMesh::Number>,QUESO::GslSparseMatrix<libMesh::Number>>
    qoiSpace(env, "qoi_", 1, NULL);

  //step 4: instantiate the qoi function object
  qoiRoutine_DataType qoiRoutine_Data;
  qoiRoutine_Data.coef1 = 1.;
  QUESO::GenericVectorFunction<QUESO::GslNumericVector<libMesh::Number>,QUESO::GslSparseMatrix<libMesh::Number>,
                               QUESO::GslNumericVector<libMesh::Number>,QUESO::GslSparseMatrix<libMesh::Number>>
    qoiFunctionObj("qoi_",
        paramDomain,
        qoiSpace,
        qoiRoutine,
        (void *) &qoiRoutine_Data);

  //step 5: instantiate the forward problem
  //parameter is Jeffreys RV
  QUESO::JeffreysVectorRV<QUESO::GslNumericVector<libMesh::Number>,QUESO::GslSparseMatrix<libMesh::Number>>
    paramRv("param_", paramDomain);

  QUESO::GenericVectorRV<QUESO::GslNumericVector<libMesh::Number>,QUESO::GslSparseMatrix<libMesh::Number>>
    qoiRv("qoi_",qoiSpace);

  QUESO::StatisticalForwardProblem<QUESO::GslNumericVector<libMesh::Number>,QUESO::GslSparseMatrix<libMesh::Number>,
                                   QUESO::GslNumericVector<libMesh::Number>,QUESO::GslSparseMatrix<libMesh::Number>>
    fp("", NULL, paramRv, qoiFunctionObj, qoiRv);

  //step 6: solve the forward problem
  fp.solveWithMonteCarlo(NULL);
}

int main(int argc, char ** argv) {
  std::string inputFileName = "test_Regression/jeffreys_input.txt";
  const char * test_srcdir = std::getenv("srcdir");
  if (test_srcdir)
    inputFileName = test_srcdir + ('/' + inputFileName);

  //init env
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc,&argv);
  // Step 1: Set up QUESO environment
  QUESO::FullEnvironment* env =
    new QUESO::FullEnvironment(MPI_COMM_WORLD,inputFileName,"",NULL);
#else
  QUESO::FullEnvironment* env =
    new QUESO::FullEnvironment(inputFileName,"",NULL);
#endif

  //compute
  compute(*env);

  //finalize enviroment
  delete env;
#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
