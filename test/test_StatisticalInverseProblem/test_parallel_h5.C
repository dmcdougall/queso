#include <queso/Environment.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>
#include <queso/GslNumericVector.h>
#include <queso/GslSparseMatrix.h>

template<class V = QUESO::GslNumericVector<libMesh::Number>, class M = QUESO::GslSparseMatrix<libMesh::Number>>
class Likelihood : public QUESO::BaseScalarFunction<V, M>
{
public:

  Likelihood(const char * prefix, const QUESO::VectorSet<V, M> & domain)
    : QUESO::BaseScalarFunction<V, M>(prefix, domain)
  {
  }

  virtual ~Likelihood()
  {
  }

  virtual double lnValue(const V & /* domainVector */, const V * /* domainDirection */,
      V * /* gradVector */, M * /* hessianMatrix */, V * /* hessianEffect */) const
  {
    double misfit = 1.0;

    return -0.5 * misfit;
  }

  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    return std::exp(this->lnValue(domainVector, domainDirection, gradVector,
          hessianMatrix, hessianEffect));
  }

  using QUESO::BaseScalarFunction<V, M>::lnValue;
};

int main(int argc, char ** argv) {
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
  QUESO::FullEnvironment env(MPI_COMM_WORLD, argv[1], "", NULL);
#else
  QUESO::FullEnvironment env(argv[1], "", NULL);
#endif

  QUESO::VectorSpace<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > paramSpace(env, "param_", 1, NULL);

  QUESO::GslNumericVector<libMesh::Number> paramMins(paramSpace.zeroVector());
  QUESO::GslNumericVector<libMesh::Number> paramMaxs(paramSpace.zeroVector());

  double min_val = 0.0;
  double max_val = 1.0;
  paramMins.cwSet(min_val);
  paramMaxs.cwSet(max_val);

  QUESO::BoxSubset<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > paramDomain("param_", paramSpace, paramMins, paramMaxs);

  QUESO::UniformVectorRV<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > priorRv("prior_", paramDomain);

  Likelihood<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > lhood("llhd_", paramDomain);

  QUESO::GenericVectorRV<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > postRv("post_", paramSpace);

  QUESO::StatisticalInverseProblem<QUESO::GslNumericVector<libMesh::Number>, QUESO::GslSparseMatrix<libMesh::Number> > ip("", NULL, priorRv, lhood, postRv);

  QUESO::GslNumericVector<libMesh::Number> paramInitials(paramSpace.zeroVector());
  paramInitials[0] = 0.0;

  QUESO::GslSparseMatrix<libMesh::Number> proposalCovMatrix(paramSpace.zeroVector());
  proposalCovMatrix(0, 0) = 0.1;

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
