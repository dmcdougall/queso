//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_SCALAR_FUNCTION_H__
#define __UQ_SCALAR_FUNCTION_H__

#include <queso/VectorSet.h>
#include <queso/VectorSubset.h>
#include <queso/Environment.h>
#include <queso/Defines.h>

namespace QUESO {

//*****************************************************
// Base class
//*****************************************************

/*! \file uqScalarFunction.h
 * \brief Set of classes for handling vector functions.
 * 
 * \class BaseScalarFunction
 * \brief A templated (base) class for handling scalar functions.
 *
 * This class allows the mathematical definition of a scalar function such as:
 * \f$ f: B \subset R \rightarrow R \f$. A function of one or more variables 
 * has always one-dimensional range.  PDFs (marginal, joint) and CDFs are examples 
 * of scalar functions. */


template<class V,class M>
class BaseScalarFunction {
public:
  //! @name Constructor/Destructor methods.
  //@{ 
  //! Default constructor.
  /*! Instantiates an object of the class, i.e. a scalar function, given a prefix and its domain.*/
  BaseScalarFunction(const char*                  prefix,
			    const VectorSet<V,M>& domainSet);
  //! Destructor	
  virtual ~BaseScalarFunction();
  //@}
  
    //! @name Mathematical methods.
  //@{  
  //! Access to the protected attribute \c m_domainSet: domain set of the scalar function.
  const VectorSet<V,M>& domainSet  ()  const;
  
  
  //! Actual value of the scalar function.
  virtual       double                 actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const = 0;
  
  //! Logarithm of the value of the scalar function.
  virtual       double                 lnValue    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const = 0;
  //@}
protected:
  const BaseEnvironment& m_env;
        std::string             m_prefix;
	
  //! Domain set of the scalar function.	
  const VectorSet<V,M>&  m_domainSet;
};
// Default constructor -----------------------------
template<class V,class M>
BaseScalarFunction<V,M>::BaseScalarFunction(
  const char*                  prefix,
  const VectorSet<V,M>& domainSet)
  :
  m_env      (domainSet.env()),
  m_prefix   ((std::string)(prefix)+"func_"),
  m_domainSet(domainSet)
{
}
// Destructor ---------------------------------------
template<class V,class M>
BaseScalarFunction<V,M>::~BaseScalarFunction()
{
}
// Math methods -------------------------------------
template<class V,class M>
const VectorSet<V,M>&
BaseScalarFunction<V,M>::domainSet() const
{
  return m_domainSet;
}

//*****************************************************
// Generic class
//*****************************************************

/*!\class GenericScalarFunction
 * \brief A class for handling generic scalar functions.
 *
 * This class allows the mathematical definition of a scalar function such as:
 * \f$ f: B \subset R \rightarrow R \f$. It is derived from BaseScalarFunction. */

template<class V,class M>
class GenericScalarFunction : public BaseScalarFunction<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor. 
  /*! Instantiates an object of \c this class given its prefix, domain set and a pointer to a routine. 
   This routine plays the role of a scalar math function, and it is useful, for instance, to calculate 
   the likelihood (and its image set).*/
  GenericScalarFunction(const char*                  prefix,
                               const VectorSet<V,M>& domainSet,
                               double (*valueRoutinePtr)(const V& domainVector, const V* domainDirection, const void* routinesDataPtr, V* gradVector, M* hessianMatrix, V* hessianEffect),
                               const void* routinesDataPtr,
                               bool routineIsForLn);
  //! Virtual destructor
  virtual ~GenericScalarFunction();

  //! @name Mathematical method
  //@{ 
  //! Calculates the actual value of this scalar function.
  double actualValue      (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  
  //! Calculates the logarithm of value of this scalar function.
  /*! It is used in routines that calculate the likelihood and expect the logarithm of value.*/
  double lnValue          (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  //@}
protected:
  using BaseScalarFunction<V,M>::m_env;
  using BaseScalarFunction<V,M>::m_prefix;
  using BaseScalarFunction<V,M>::m_domainSet;

   //! Routine defining a scalar function.
   /*! The presence of the parameters \c gradVectors, \c hessianMatrices and \c hessianEffects
   * allows the user to calculate gradient vectors, Hessian matrices and Hessian effects; which 
   * can hold important information about her/his statistical application. Used, for instance to 
   * define the likelihood.  */
  double (*m_valueRoutinePtr)(const V& domainVector, const V* domainDirection, const void* routinesDataPtr, V* gradVector, M* hessianMatrix, V* hessianEffect);
  const void* m_routinesDataPtr;
  bool m_routineIsForLn;
};
// Default constructor -----------------------------
template<class V,class M>
GenericScalarFunction<V,M>::GenericScalarFunction(
  const char*                  prefix,
  const VectorSet<V,M>& domainSet,
  double (*valueRoutinePtr)(const V& domainVector, const V* domainDirection, const void* routinesDataPtr, V* gradVector, M* hessianMatrix, V* hessianEffect),
  const void* routinesDataPtr,
  bool routineIsForLn)
  :
  BaseScalarFunction<V,M>(((std::string)(prefix)+"gen").c_str(), domainSet),
  m_valueRoutinePtr             (valueRoutinePtr),
  m_routinesDataPtr             (routinesDataPtr),
  m_routineIsForLn              (routineIsForLn)
{
}
// Destructor ---------------------------------------
template<class V,class M>
GenericScalarFunction<V,M>::~GenericScalarFunction()
{
}
// Math methods -------------------------------------
template<class V,class M>
double
GenericScalarFunction<V,M>::actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO(m_valueRoutinePtr == NULL,
                      m_env.worldRank(),
                      "GenericScalarFunction<V,M>::actualValue()",
                      "m_valueRoutinePtr = NULL");

  double value = m_valueRoutinePtr(domainVector, domainDirection, m_routinesDataPtr, gradVector, hessianMatrix, hessianEffect);
  if (m_routineIsForLn) {
#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
    value = std::exp(value);
#else
    value = std::exp(-.5*value);
#endif
    UQ_FATAL_TEST_MACRO((domainDirection != NULL) ||
                        (gradVector      != NULL) ||
                        (hessianMatrix   != NULL) ||
                        (hessianEffect   != NULL),
                        m_env.worldRank(),
                        "GenericScalarFunction<V,M>::gradOfActual()",
                        "INCOMPLETE CODE");
  }
  return value;
}
// --------------------------------------------------
template<class V,class M>
double
GenericScalarFunction<V,M>::lnValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO(m_valueRoutinePtr == NULL,
                      m_env.worldRank(),
                      "GenericScalarFunction<V,M>::lnValue()",
                      "m_valueRoutinePtr = NULL");

  double value = m_valueRoutinePtr(domainVector, domainDirection, m_routinesDataPtr, gradVector, hessianMatrix, hessianEffect);
  if (m_routineIsForLn == false) {
#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
    value = log(value);
#else
    value = -2.*log(value);
#endif
    UQ_FATAL_TEST_MACRO((domainDirection != NULL) ||
                        (gradVector      != NULL) ||
                        (hessianMatrix   != NULL) ||
                        (hessianEffect   != NULL),
                        m_env.worldRank(),
                        "GenericScalarFunction<V,M>::gradOfLn()",
                        "INCOMPLETE CODE");
  }
  return value;
}

//*****************************************************
// Constant class
//*****************************************************

/*!\class ConstantScalarFunction
 * \brief A class for handling scalar functions which image is a constant (real number).
 *
 * This class allows the mathematical definition of a scalar function which image set 
 * is a constant (real number). */

template<class V,class M>
class ConstantScalarFunction : public BaseScalarFunction<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! Instantiates an object of the class, i.e. a scalar function, given a prefix, its domain and constant-valued image.*/
  ConstantScalarFunction(const char*                  prefix,
                                const VectorSet<V,M>& domainSet,
                                double                       constantValue);
  //! Virtual destructor
  virtual ~ConstantScalarFunction();

  //! @name Mathematical method
  //@{ 
  //! Calculates the actual value of this scalar function.
  double actualValue      (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  
  //! Calculates the logarithm of the value of this scalar function (which is zero).
  double lnValue          (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  //@}
protected:
  using BaseScalarFunction<V,M>::m_env;
  using BaseScalarFunction<V,M>::m_prefix;
  using BaseScalarFunction<V,M>::m_domainSet;

  //! Constant value is the image set of this scalar function.
  double m_constantValue;
};
// Default constructor -----------------------------
template<class V,class M>
ConstantScalarFunction<V,M>::ConstantScalarFunction(
  const char*                  prefix,
  const VectorSet<V,M>& domainSet,
  double                       constantValue)
  :
  BaseScalarFunction<V,M>(((std::string)(prefix)+"gen").c_str(), domainSet),
  m_constantValue               (constantValue)
{
}
// Destructor ---------------------------------------
template<class V,class M>
ConstantScalarFunction<V,M>::~ConstantScalarFunction()
{
}
// Math methods -------------------------------------
template<class V,class M>
double
ConstantScalarFunction<V,M>::actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const
{
  return m_constantValue;
}
// --------------------------------------------------
template<class V,class M>
double
ConstantScalarFunction<V,M>::lnValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const
{
  return 0.;
}

}  // End namespace QUESO

#endif // __UQ_SCALAR_FUNCTION_H__
