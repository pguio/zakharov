/**************************************************************************
 *
 * $Id: mapper.h,v 1.68 2014/10/21 10:25:44 patrick Exp $
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

#ifndef MAPPER_H
#define MAPPER_H

#include <zakharov-defs.h>
#include <SIconsts.h>
#include <fourier.h>

class Mapper : public parser::Parser {

  typedef std::ostream ostream;
  typedef std::vector<Real> Vector;

  typedef blitz::Range Range;

  enum parser_enum {__N=1, __M,
                    _physical_units, __L, __h,
                    _ne, _Te, _nuec, _Ne,
                    _nb_ne, _ub, _dub_ub,
                    _mi, _Ti, _nuic,
                    _dOmega, _Epump
                   };
public:

  friend ostream& operator<<(ostream& os, const Mapper &mapper);

  Mapper(int nargs, char *args[]);
  ~Mapper();

  void initialise();
  int getGridSize() const {
    return N;
  }
  int getZeroPaddingSize() const {
    return M;
  }
  Real getPhysicalTimeStep() const {
    return toPhysicalTime(h);
  }
  Real getNormalisedTimeStep() const {
    return h;
  }
  Real getPhysicalSystemLength() const {
    return toPhysicalLength(L);
  }
  Real getNormalisedSystemLength() const {
    return L;
  }
  Real getDensity() const {
    return ne;
  }
  Array1dr getPhysicalWaveNumbers() const;
  Array1dr getNormalisedWaveNumbers() const;
  Array1dr getAlphas() const;
  Array1dr getPhysicalGridPoints() const;
  Array1dr getNormalisedGridPoints() const;
  Array1dr getDiffusionWavenumbers() const;
  Array1dr getDiffusionAlphas() const;
  Array1dr getDiffusionElectronPdf() const;
  Array1dr getDiffusionLangmuirWavesDamping() const;
  Array1dr getnkFluctuationsLevel() const;
  Array1dr getEkFluctuationsLevel() const;
  Array1dr getSoundWavesDamping() const;
  Array1dr getMaxwellianLangmuirWavesDamping() const;
  Array1dr getLangmuirWavesDamping() const;

  Real getPhysicalElectronCollisionFreq() const {
    return nuec;
  }
  Real getNormalisedElectronCollisionFreq() const {
    return toNormalisedFreq(nuec);
  }

  Complex getPump(Real t) const;
  Real getN2Wconstant() const;
  Complex getDispersionRelation(Real k, Complex w) const;
  Complex getDispersionRelation1stDerivative(Real k, Complex w) const;
  Complex getDispersionRelation2ndDerivative(Real k, Complex w) const;
  Complex getDispersionRelationRoot(Real k, Complex wi) const;
  Complex getDispersionRelationRoot1(Real k, Complex wi) const;
  Array1dr initDiffusion();
  void updateDiffusion(const Array1dc &Ek, Array1dr &nue);

  template<class T>
  T toPhysicalTime(T t) const {
    return T(t*tau);
  }
  template<class T>
  T toPhysicalFreq(T f) const {
    return T(f*(1.0/tau));
  }
  template<class T>
  T toPhysicalLength(T l) const {
    return T(l*chi);
  }
  template<class T>
  T toPhysicalWaveNumber(T k) const {
    return T(k*(1.0/chi));
  }
  template<class T>
  T toPhysicalElectricField(T E) const {
    return T(E*this->epsilon);
  }
  template<class T>
  T toPhysicalDensity(T n) const {
    return T(n*nu);
  }
  template<class T>
  T toPhysicalWaveNumber2(T k) const {
    return T(k*(1.0/blitz::pow2(chi)));
  }
  template<class T>
  T toPhysicalElectricField2(T E2) const {
    return T(E2*blitz::pow2(this->epsilon));
  }
  template<class T>
  T toPhysicalDensity2(T n2) const {
    return T(n2*blitz::pow2(nu));
  }
  template<class T>
  T toPhysicalPlasmaLine(T s) const {
#if 0
    return T(s*(SI::me/mi)*blitz::pow2(SI::eps0/(SI::Kb*Te)));
#endif

    return T(s*blitz::pow2(SI::eps0/SI::e));
  }

  template<class T>
  T toNormalisedTime(T t) const {
    return T(t*(1.0/tau));
  }
  template<class T>
  T toNormalisedFreq(T f) const {
    return T(f*tau);
  }
  template<class T>
  T toNormalisedLength(T l) const {
    return T(l*(1.0/chi));
  }
  template<class T>
  T toNormalisedWaveNumber(T k) const {
    return T(k*chi);
  }
  template<class T>
  T toNormalisedElectricField(T E) const {
    return T(E*(1.0/this->epsilon));
  }
  template<class T>
  T toNormalisedDensity(T n) const {
    return T(n*(1.0/nu));
  }

protected:

  typedef parser::Parser Parser;

private:
  int N;
  int M;

  bool physical_units;

  // Physical parameters
  Real L_m;      // [m]
  Real h_s;      // [s]

  // Normalised parameters
  Real L;        // []
  Real h;        // []

  // Physical input parameters
  Real ne;       // [m-3]
  Real Te;       // [K]
  Real nuec;     // [s-1]
  Real Ne;       // []

  Vector nb_ne;   // []
  Vector ub;      // [m s-1]
  Vector dub_ub;  // []
  Vector Eb;      // [eV]
  Vector kb;      // [m-1]
  Vector Gam;     // \Gamma = (nb/ne)*(vb/delta vb)
  Vector W_L;     // time-asymptotic energy in the electron beam instability
  Vector W_Lc;    // time-asymptotic energy in the cold electron beam instability

  Real mni;      // [amu]
  Real mi;       // [kg]
  Real Ti;       // [K]
  Real nuic;     // [s-1]

  // Physical pump parameters
  Real dOmega;   // [rad-1] dOmega= wpump - wpe
  Real Epump;    // [V m-1]
  Real wpump;    // [rad-1]
  Real kpump;    // [m-1]

  // Normalised pump parameters
  Real dOmega_n; // []
  Real Epump_n;  // []
  Real kpump_n;  // []

  // Plasma parameters
  Real eta;      // [] eta = (ce*Te+ci*Ti)/Te
  Real wpe;      // [rad s-1]
  Real ve;       // [m s-1]
  Real vi;       // [m s-1]
  Real cs;       // [m s-1]
  Real ve2;      // [m s-1]^2
  Real lambdae;  // [m]
  Real ke;       // [m-1]

  // Zakharov equation normalisation constants
  Real tau;      // [s]
  Real chi;      // [m]
  Real epsilon;  // [V m-1]
  Real nu;       // [m-3]

  // Diffusion variable
  Array1dr vphi;
  Array1dr Fet0;
//#define SEP 1
#if SEP // Fem separated
  Array1dr Fem;
#endif
  Array1dr Fe;
  Array1dr dFe;
  Array1dr dkdvphi;
  Array1dr dampingCoeff;
  Array1dr gamma;
  Array1dr D;
  Array1dr dD;
  Array1dr diffusionCoeff;
#if 0 // alternative diffusion term

  Array1dr wk, gk, gk2;
#endif
#if 0 // buffer to estimate time mean value of <E(k)^2>

  Array2dr Ekt2;
#endif

  Range ipRange;
  Range inRange;
  Range ip1Range;
  Range in1Range;
  Real dk;
#if 0 // buffer to estimate time mean value of <E(k)^2>

  int nt;
  int it;
  bool bufferComplete;
#endif

  void initParsing(int nargs, char *args[]);
  void paramParsing();
  void checkParam() const;

  void setIntegrationBoundaries();
  void deriveWithK(const Array1dr &y, Array1dr &dy);
  void deriveFeWithVphi(Range &range);
  void ftcsLax(Array1dr &u, Array1dr &d, Array1dr &dkdv);
  void ftcsLax(Range &range);
  void crankNicholson(Array1dr &u, Array1dr &d, Array1dr &dkdv);
  void crankNicholson(Array1dr &u, Array1dr &d, Array1dr &dkdv,
                      Array1dr &v, Array1dr &fet0);
  void crankNicholson(Range &range);
  void tridiag(Array1dr &a, Array1dr &b, Array1dr &c, Array1dr &r, Array1dr &u);
  Complex wofz(const Complex z) const;
  Complex W(const Complex z) const;
  Complex dWdz(const Complex z) const;
  Complex d2Wdz2(const Complex z) const;
#if defined(HAVE_GNUPLOT)
  void gnuplot();
#endif
};


#endif // MAPPER_H

