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

#include <queso/Fft.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_complex.h>

#ifdef QUESO_HAVE_EIGEN
#include <unsupported/Eigen/FFT>
#endif  // QUESO_HAVE_EIGEN

namespace QUESO {

#ifdef QUESO_HAVE_EIGEN
void
eigen_impl_fft_real_forward(const std::vector<double> & data,
                       unsigned int fftSize,
                       std::vector<std::complex<double> > & forwardResult)
{
  Eigen::FFT<double> fft;
  fft.fwd(forwardResult, data);
}
#else
void
gsl_impl_fft_real_forward(const std::vector<double> & data,
                     unsigned int fftSize,
                     std::vector<std::complex<double> > & forwardResult)
{
  gsl_fft_real_workspace* realWkSpace = gsl_fft_real_workspace_alloc(fftSize);
  gsl_fft_real_wavetable* realWvTable = gsl_fft_real_wavetable_alloc(fftSize);

  std::vector<double> data_tmp(data);
  gsl_fft_real_transform(&data_tmp[0],
                         1,
                         fftSize,
                         realWvTable,
                         realWkSpace);

  gsl_fft_real_wavetable_free(realWvTable);
  gsl_fft_real_workspace_free(realWkSpace);

  unsigned int halfFFTSize = fftSize/2;
  bool sizeIsEven = ((fftSize % 2) == 0);
  double realPartOfFFT = 0.;
  double imagPartOfFFT = 0.;
  for (unsigned int j = 0; j < data_tmp.size(); ++j) {
    if (j == 0) {
      realPartOfFFT = data_tmp[j];
      imagPartOfFFT = 0.;
    }
    else if (j < halfFFTSize) {
      realPartOfFFT = data_tmp[2*j-1];
      imagPartOfFFT = data_tmp[2*j  ];
    }
    else if (j == halfFFTSize) {
      realPartOfFFT = data_tmp[2*j-1];
      if (sizeIsEven) {
        imagPartOfFFT = 0.;
      }
      else {
        imagPartOfFFT = data_tmp[2*j  ];
      }
    }
    else {
      realPartOfFFT =  data_tmp[2*(fftSize-j)-1];
      imagPartOfFFT = -data_tmp[2*(fftSize-j)  ];
    }
    forwardResult[j] = std::complex<double>(realPartOfFFT,imagPartOfFFT);
  }
}
#endif  // QUESO_HAVE_EIGEN

void
eigen_impl_fft_real_inverse(const std::vector<double> & data,
                            unsigned int fftSize,
                            std::vector<std::complex<double> > & inverseResult)
{
  std::vector<std::complex<double> > internalData(fftSize,0.);
  unsigned int minSize = std::min((unsigned int) data.size(),fftSize);
  for (unsigned int j = 0; j < minSize; ++j) {
    internalData[j] = std::complex<double>(data[j], 0.0);
  }

  Eigen::FFT<double> fft;
  fft.inv(inverseResult, internalData);
}

void
gsl_impl_fft_real_inverse(const std::vector<double> & data,
                          unsigned int fftSize,
                          std::vector<std::complex<double> > & inverseResult)
{
  std::vector<double> internalData(2*fftSize,0.); // Yes, twice the fftSize
  unsigned int minSize = std::min((unsigned int) data.size(),fftSize);
  for (unsigned int j = 0; j < minSize; ++j) {
    internalData[2*j] = data[j];
  }

  gsl_fft_complex_workspace* complexWkSpace = gsl_fft_complex_workspace_alloc(fftSize);
  gsl_fft_complex_wavetable* complexWvTable = gsl_fft_complex_wavetable_alloc(fftSize);

  gsl_fft_complex_inverse(&internalData[0],
                          1,
                          fftSize,
                          complexWvTable,
                          complexWkSpace);

  gsl_fft_complex_wavetable_free(complexWvTable);
  gsl_fft_complex_workspace_free(complexWkSpace);

  for (unsigned int j = 0; j < fftSize; ++j) {
    inverseResult[j] = std::complex<double>(internalData[2*j],internalData[2*j+1]);
  }
}

// Math methods------------------------------------------
template <>
void
Fft<double>::forward(
  const std::vector<double>&                data,
        unsigned int                        fftSize,
        std::vector<std::complex<double> >& forwardResult)
{
  if (forwardResult.size() != fftSize) {
    forwardResult.resize(fftSize,std::complex<double>(0.,0.));
    std::vector<std::complex<double> >(forwardResult).swap(forwardResult);
  }

  std::vector<double> internalData(fftSize,0.);
  unsigned int minSize = std::min((unsigned int) data.size(),fftSize);
  for (unsigned int j = 0; j < minSize; ++j) {
    internalData[j] = data[j];
  }

#ifdef QUESO_HAVE_EIGEN
  eigen_impl_fft_real_forward(internalData, fftSize, forwardResult);
#else
  gsl_impl_fft_real_forward(internalData, fftSize, forwardResult);
#endif
}
//-------------------------------------------------------
template <>
void
Fft<double>::inverse(
  const std::vector<double>&                data,
        unsigned int                        fftSize,
        std::vector<std::complex<double> >& inverseResult)
{
  if (inverseResult.size() != fftSize) {
    inverseResult.resize(fftSize,std::complex<double>(0.,0.));
    std::vector<std::complex<double> >(inverseResult).swap(inverseResult);
  }

  gsl_impl_fft_real_inverse(data, fftSize, inverseResult);
}

}  // End namespace QUESO
