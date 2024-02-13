/**************************************************************************
 *
 * $Id: zakharov-defs.h,v 1.56 2019/05/10 16:54:43 patrick Exp $
 *
 * Copyright (c) 2001-2011 Patrick Guio <patrick.guio@gmail.com>
 * All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2.  of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef ZAKHAROV_DEFS_H
#define ZAKHAROV_DEFS_H

#if defined(HAVE_CONFIG_H)
#include <zakharov-config.h>
#endif

#include <list>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>

#include <blitz/array.h>
#include <blitz/numinquire.h>
#include <random/uniform.h>

#include "parser.h"

#include "range-spec.h"

#if defined(HAVE_PLPLOT)
#include <plstream.h>
#endif

#define ZAK_COPYRIGHT \
"Copyright (c) 2001-2011 Patrick Guio <patrick.guio@gmail.com>\n\n"\
"This is free software; see the source for copying conditions.  There is NO\n"\
"warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"

typedef double Real;
typedef std::complex<Real> Complex;

typedef blitz::Array<Complex, 1> Array1dc;
typedef blitz::Array<Complex, 2> Array2dc;
typedef blitz::Array<Complex, 3> Array3dc;

typedef blitz::Array<Real, 1> Array1dr;
typedef blitz::Array<Real, 2> Array2dr;
typedef blitz::Array<Real, 3> Array3dr;

typedef blitz::Array<int, 1> Array1di;

typedef blitz::TinyVector<int,2> Vector2di;
typedef blitz::TinyVector<int,3> Vector3di;

const RangeSpec<int> DEFAULT_SaveIter(0, 10, 50);

template<class T>
std::ostream& operator<<(std::ostream &os, const std::vector<T> &v)
{
  if (!v.empty()) {
    typename std::vector<T>::const_iterator i = v.begin();
    typename std::vector<T>::const_iterator end = v.end();
    os << *i;
    for ( ++i; i != end; ++i ) os << ',' << *i;
  }
  return os;
}

// Miscelleneous constant values
const Real m_pi      = 3.14159265358979323846;  // pi
const Real m_2pi     = 2.0*m_pi;                // 2 pi
const Real m_sqrt1_2 = 0.70710678118654752440;  // 1/sqrt(2)
const Real m_sqrt_pi = 1.77245385090551603;     // sqrt(pi)
const Complex I      = Complex(0,1);            // sqrt(-1)

// Default solver parameters
const int NFFT_DEFAULT = 64;
const Real L_DEFAULT   = 20.0;
const Real H_DEFAULT   = 5e-3;

inline
Real integrate1d(Array1dr &x, Array1dr &f)
{
  int n(x.size());
  blitz::Range i(0,n-2);
#if 1
  // compount trapezoidal method
  Array1dr fdx(0.5*(f(i)+f(i+1))*(x(i+1)-x(i)));
#else
  Array1dr fdx(f*(x(1)-x(0)));
#endif
  return blitz::sum(fdx);
}

inline
Array1dc complexExp(const Array1dc & z)
{
  using blitz::exp;
#ifndef BZ_HAVE_COMPLEX_MATH
  return Array1dc(exp(real(z))*zip(cos(imag(z)), sin(imag(z)), Complex()));
#else
  return Array1dc(exp(z));
#endif
}

inline
Array1dr AbsSquare(const Array1dc & z)
{
  return Array1dr(pow2(real(z))+pow2(imag(z)));
}

template<class T>
void saveMatlab(const std::string filename, T mode,
                const std::string varname, const Real var)
{
  std::ofstream os(filename.c_str(), mode);
  os.setf(std::ios::scientific);
  os.precision(6);
  os << varname << " = " << var << ";\n";
}

template<class T>
void saveMatlab(const std::string filename, T mode,
                const std::string varname, const Array1dr var)
{
  std::ofstream os(filename.c_str(), mode);
  os.setf(std::ios::scientific);
  os.precision(6);
  os << varname << " = [ ";
  for (unsigned i=0; i<var.size(); ++i) {
    os << var(i) << "; ";
    if (i != 0 && (i%10) == 0)  os << "...\n";
  }
  os << "];\n";
}

template<class T>
void saveMatlab(const std::string filename, T mode,
                const std::string varname, const Array1dc var)
{
  std::ofstream os(filename.c_str(), mode);
  os.setf(std::ios::scientific);
  os.precision(6);
  os << varname << " = [ ";
  for (unsigned i=0; i<var.size(); ++i) {
    os << std::real(var(i))
       << "+i*" << std::imag(var(i)) << "; ";
    if (i != 0 && (i%10) == 0)  os << "...\n";
  }
  os << "];\n";
}

inline
int adjfn(int st, int fn, int stride)
{
  return fn - (fn-st) % stride;
}

inline
int checkNaN(const std::string classname, const Array1dc & var,
             const std::string varname, const std::string filename,
             int linenumber)
{
  using blitz::blitz_isnan;
  using blitz::any;
  if ( any(blitz_isnan(real(var))) || any(blitz_isnan(imag(var))) ) {
    std::ostringstream os;
    os << "NaN detected in " << varname << " in file: " << filename << ": "
       << linenumber << ".";
    throw ClassException(classname, os.str());
  }
  return 0;
}

inline
void printOsStatus(std::ostream &os)
{
  std::ios::fmtflags f = os.flags();
  int p = os.precision();
  int w = os.width();
  char c = os.fill();
  std::cout << "f=" << f << " p=" << p << " c='" << c << "' w=" << w << std::endl;
}

#endif // ZAKHAROV_DEFS_H

