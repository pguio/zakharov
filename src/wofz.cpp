
#include <iostream>
#include <cmath>
#include <complex>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

#define SQRT_PI 1.77245385090551603

typedef std::complex<double> Complex;


#define I Complex(0,1)

template<typename Real>
std::complex<double> wofz(const std::complex<Real> &z)
{
  // Algorithm 680, collected algorithms from acm.
  // This work published in transactions on mathematical software,
  // vol. 16, no. 1, pp. 47.


  // Given a complex number z = (xi,yi), this subroutine computes
  // the value of the Faddeeva-function w(z) = exp(-z^2)*erfc(-iz),
  // where erfc is the complex complementary error-function and i
  // means sqrt(-1).
  // The accuracy of the algorithm for z in the 1st and 2nd quadrant
  // is 14 significant digits; in the 3rd and 4th it is 13 significant
  // digits outside a circular region with radius 0.126 around a zero
  // of the function.
  // All real variables in the program are double precision.

  // The code contains a few compiler-dependent parameters :
  //     rmaxreal = the maximum value of rmaxreal equals the root oF
  //                rmax = the largest number which can still be
  //                implemented on the computer in double precision
  //                floating-point arithmetic
  //     rmaxexp  = ln(rmax) - ln(2)
  //     rmaxgoni = the largest possible argument of a double precision
  //                goniometric function (cos, sin, ...)
  // The reason why these parameters are needed as they are defined will
  // be explained in the code by means of comments
  //
  // Reference - GPM Poppe, CMJ Wijers; more efficient computation oF
  // the complex error-function, acm trans. math. software.
  Real xi = real(z);
  Real yi = imag(z);
  Real u;
  Real v;

  const Real factor = 1.12837916709551257388; // 2/sqrt(pi)
  const Real rmaxreal = 0.5e154;
  const Real rmaxexp = 708.503061461606;
  const Real rmaxgoni = 3.53711887601422e15;

  Real xabs = fabs(xi);
  Real yabs = fabs(yi);
  Real x = xabs/6.3;
  Real y = yabs/4.4;

  // The following if-statement protects
  // qrho = (x^2 + y^2) against overflow
  if (xabs>rmaxreal || yabs>rmaxreal) {
    std::cerr << "Error in wofz\n\txabs>rmaxreal || yabs>rmaxreal\n";
    throw !0;
  }

  Real qrho = x*x + y*y;

  Real xabsq = xabs*xabs;
  Real xquad = xabsq - yabs*yabs;
  Real yquad = 2.0*xabs*yabs;

  Real u2;
  Real v2;

  bool a = (qrho < 0.085264);

  if (a) {
    // If (qrho < 0.085264) then the Faddeeva-function is evaluated
    // using a power-series (Abramowitz/Stegun, equation (7.1.5), p.297)
    // n is the minimum number of terms needed to obtain the required
    // accuracy
    qrho  = (1.0-0.85*y)*sqrt(qrho);
    int n = int(6.0 + 72.0*qrho);
    int j = 2*n+1;
    Real xsum  = 1.0/Real(j);
    Real ysum  = 0.0;
    for (int i=n; i>=1; --i) {
      j = j - 2;
      Real xaux = (xsum*xquad - ysum*yquad)/Real(i);
      ysum = (xsum*yquad + ysum*xquad)/Real(i);
      xsum = xaux + 1.0/Real(j);
    }
    Real u1   = -factor*(xsum*yabs + ysum*xabs) + 1.0;
    Real v1   =  factor*(xsum*xabs - ysum*yabs);
    Real daux = exp(-xquad);
    u2   =  daux*cos(yquad);
    v2   = -daux*sin(yquad);

    u = u1*u2 - v1*v2;
    v = u1*v2 + v1*u2;

  } else {
    // If (qrho>1.0) then w(z) is evaluated using the Laplace
    // continued fraction
    // nu is the minimum number of terms needed to obtain the required
    // accuracy

    // If (qrho > 0.085264) && qrho < 1.0) then w(z) is evaluated
    // by a truncated Taylor expansion, where the Laplace continued fraction
    // is used to calculate the derivatives of w(z)
    // kapn is the minimum number of terms in the taylor expansion needed
    // to obtain the required accuracy
    // nu is the minimum number of terms of the continued fraction needed
    // to calculate the derivatives with the required accuracy
    Real h;
    Real h2;
    int kapn;
    int nu;
    if (qrho > 1.0) {
      h = 0.0;
      kapn = 0;
      qrho = sqrt(qrho);
      nu = int(3.0 + (1442.0/(26.0*qrho+77.0)));
    } else {
      qrho = (1-y)*sqrt(1.0-qrho);
      h = 1.88*qrho;
      h2 = 2.0*h;
      kapn = int(7.0 + 34.0*qrho);
      nu = int(16.0 + 26.0*qrho);
    }

    bool b = (h > 0.0);

    Real qlambda;
    if (b)
      qlambda = pow(h2, Real(kapn));

    Real rx = 0.0;
    Real ry = 0.0;
    Real sx = 0.0;
    Real sy = 0.0;

    for (int n=nu; n>=0; --n) {
      int np1 = n + 1;
      Real tx  = yabs + h + np1*rx;
      Real ty  = xabs - np1*ry;
      Real c   = 0.5/(tx*tx + ty*ty);
      rx  = c*tx;
      ry  = c*ty;
      if (b && n <= kapn) {
        tx = qlambda + sx;
        sx = rx*tx - ry*sy;
        sy = ry*tx + rx*sy;
        qlambda = qlambda/h2;
      }
    }

    if (h == 0.0) {
      u = factor*rx;
      v = factor*ry;
    } else {
      u = factor*sx;
      v = factor*sy;
    }

    if (yabs == 0.0)
      u = exp(-xabs*xabs);
  }

  // Evaluation of w(z) in the other quadrants
  if (yi < 0.0) {
    if (a) {
      u2    = 2*u2;
      v2    = 2*v2;
    } else {
      xquad =  -xquad;

      // The following if-statement protects 2*exp(-z^2)
      // against overflow
      if (yquad>rmaxgoni || xquad>rmaxexp) {
        std::cerr << "Error in wofz\n\txabs>rmaxreal || yabs>rmaxreal\n";
        throw !0;
      }

      Real w1 =  2.0*exp(xquad);
      u2  =  w1*cos(yquad);
      v2  = -w1*sin(yquad);
    }

    u = u2 - u;
    v = v2 - v;

    if (xi > 0.0)
      v = -v;
  } else {
    if (xi < 0.0)
      v = -v;
  }

  return std::complex<double>(u,v);
}

template<typename Real>
std::complex<Real> W(const std::complex<Real> &z)
{
  std::complex<Real> z2 = std::complex<Real>(M_SQRT1_2,0.0)*z;
  return std::complex<Real>(1.0,0.0)+std::complex<Real>(0.0,SQRT_PI)*z2*wofz(z2);
}


int main()
{
  try {
    std::complex<double> z;
    std::cout.precision(16);
    while ( std::cin >> z.real() >> z.imag()) {
      std::cout << "z=" <<  z << " w(z)=" << wofz(z) << std::endl;
      std::cout << "z=" <<  z << " W(z)=" << W(z) << std::endl;
    }
    return 0;
  } catch (int & error) {
    return error;
  }
}
