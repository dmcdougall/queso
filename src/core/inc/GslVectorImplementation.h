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
//--------------------------------------------------------------------------

#ifndef UQ_GSL_VECTOR_IMPLEMENTATION_H
#define UQ_GSL_VECTOR_IMPLEMENTATION_H

/*! \file GslVectorImplementation.h
    \brief Vector class using GSL
*/

#include <queso/Defines.h>
#include <queso/Vector.h>
#include <gsl/gsl_vector.h>

namespace QUESO {

/*! \class GslVectorImplementation

    \brief Class for vector operations using GSL library.
    
    This class creates and provides basic support for vectors of templated 
    type as a specialization of Vector using GSL vectors, which are defined 
    by an encapsulated gsl_vector structure.
*/

class GslVectorImplementation
{
public:
  
  //! @name Constructor/Destructor methods.
  //@{ 

  GslVectorImplementation(const Map& map);
  GslVectorImplementation(const Map& map, double value);
  GslVectorImplementation(double d1, double d2, const Map& map);

  //! Construct a vector, with length the same as \c v, with evenly spaced numbers from \c start to \c end, inclusive
  GslVectorImplementation(const GslVectorImplementation&         v, double start, double end);
  GslVectorImplementation(const GslVectorImplementation&         y);
  
  //! Destructor
  ~GslVectorImplementation();
  //@}
 
  //! @name Set methods.
  //@{ 
  //! 	Copies values from vector rhs to \c this. 
  GslVectorImplementation& operator=(const GslVectorImplementation& rhs);
  
  //! Stores in \c this the coordinate-wise multiplication of \c this and a.
  GslVectorImplementation& operator*=(double a);
  
  //! Stores in \c this the coordinate-wise division of \c this by a.
  GslVectorImplementation& operator/=(double a);
  
  //! Stores in \c this the coordinate-wise multiplication of \c this with rhs.
  GslVectorImplementation& operator*=(const GslVectorImplementation& rhs);
  
  //! Stores in \c this the coordinate-wise division of \c this by rhs.
  GslVectorImplementation& operator/=(const GslVectorImplementation& rhs);
  
   //! Stores in \c this the coordinate-wise addition of \c this and rhs.
  GslVectorImplementation& operator+=(const GslVectorImplementation& rhs);
  
  //! Stores in \c this the coordinate-wise subtraction of \c this by rhs.
  GslVectorImplementation& operator-=(const GslVectorImplementation& rhs);
  //@}
  
  //! @name Accessor methods.
  //@{
  //! Element access method (non-const).
  /*! Returns the i-th element if x[i] is specified, the expression x(i) will return the same element.*/  
            double& operator[](unsigned int i);
	    
  //! Element access method (const).
  /*! Returns the i-th element if x[i] is specified, the expression x(i) will return the same element.*/ 
	    
      const double& operator[](unsigned int i) const;
  //@}

  //! @name Attribute methods.
  //@{ 
  //! Returns the length of this vector.
  unsigned int sizeLocal        () const;

  //! Returns the global length of this vector.
  unsigned int sizeGlobal       () const;
  //@}
  
  //! @name Set methods.
  //@{ 
  //! Component-wise sets all values to \c this with value.
  void         cwSet            (double value);
  //@}

  //! This function sorts the elements of the vector \c this in ascending numerical order. 
  void         sort             ();
    
  // Necessary for GslMatrix::invertMultiply() and GslMatrix::setRow/Column
  gsl_vector*  data                          () const; 

//! @name Attribute methods.
  //@{ 
  //! Returns the maximum value in the vector \c this.
  double       getMaxValue      () const;
  
  //! Returns minimum value in the vector \c this.
  double       getMinValue      () const;
  
  //! This function returns the index of the maximum value in the vector \c this.
  int          getMaxValueIndex () const;
  
  //! This function returns the index of the minimum value in the vector \c this.
  int          getMinValueIndex () const;
  
  //! This function returns maximum value in the vector \c this and its the index.
  void         getMaxValueAndIndex( double& value, int& index );
  
  //! This function returns minimum value in the vector \c this and its the index.
  void         getMinValueAndIndex( double& value, int& index );

  //! This function returns absolute value of elements in \c this.
  GslVectorImplementation abs() const;
 //@}
private:

  //! This function copies the elements of the vector src into \c this.
  void         copy             (const GslVectorImplementation& src);

  //! GSL vector.
  gsl_vector* m_vec;
};

}  // End namespace QUESO

#endif // UQ_GSL_VECTOR_IMPLEMENTATION_H
