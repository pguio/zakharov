/**************************************************************************
 *
 * $Id: zakharov.h,v 1.158 2014/10/21 11:53:10 patrick Exp $
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

#ifndef ZAKHAROV_H
#define ZAKHAROV_H

#include <zakharov-defs.h>
#include <spectra.h>
#include <modes-factory.h>
#include <default-mode.h>
#include <diagnostics.h>
#include <mapper.h>
#include <time-utils.h>

class ZakharovSolver : public parser::Parser {
public:

  typedef std::string string;
  typedef std::ostream ostream;
  typedef std::ofstream ofstream;

  typedef fourier::ODFT1D<fourier::complex,fourier::complex> FT1d;

  typedef blitz::Range Range;

  friend ostream& operator<<(ostream& os, const ZakharovSolver &s);

  enum init_enum {
    _solitons, _efluct, _tfluct, _egauss, _densfluct, _densparab
  };

  enum intg_enum {
    _classic=0, _pump, _fixed_n, _fixed_E, _linear, _noFp,
    last_intg=_noFp+1
  };

  ZakharovSolver(int nargs, char *args[]);
  ~ZakharovSolver();

  bool parseHelp() const;
  bool parseVersion() const;
  bool parseTemplate() const;

  void initialise();
  void solve();
  void switchToThermal(std::ostream &os);

protected:

  typedef parser::Parser Parser;
  typedef Parser::LUT LUT;
  typedef Parser::LUTPair Pair;

private:

  typedef ModesFactory<InitMode> Factory;

  enum parser_enum {
    _init_mode, _intg_mode, _maxiter, _modulo_display,
    _nk_source, _Ek_source, _seed,
    _sound_damping, _langmuir_damping, _diffusion,
    _time_thermal, _w_thermal,
    _optfftw3, _plplot
  };


  int N, M;
  Real h, hOver2;
  Complex ihOver2;
  Real L;
  bool nk_source, Ek_source;
  Vector2di seed;
  bool sound_damping, langmuir_damping;
  bool diffusion;

  Real time_thermal;
  Real w_thermal;
  bool thermal;

  bool optfftw3;

  unsigned iter, maxiter;
  unsigned modulo_display;

  timeutils::Timer timer;
  FT1d fft;

  Mapper mapper;

  int init_mode;
  LUT init_map;
  InitMode *mode;

  int intg_mode;
  LUT intg_map;

  Spectra spectra;
  Diagnostics diagnostics;

  string logfile;
  ofstream logs;

#if defined(HAVE_PLPLOT)
  plstream *plsnEx, *plsnEk, *plsdiffk;
  bool plplot;
  PLFLT *xp, *yp;
  PLFLT xpmin, xpmax, ypmin, ypmax;
#endif

  // wavenumber index ranges for zero padding,
  // positive k's and negative k's
  Range k0, kp, kn;

  // Grid variables ending with
  // j : configuration space variables
  // k : spectral variables
  Array1dr xj, k, k2, hOver2k2;
  Array1dc ik2;
  Array1dc Ej, Ek;
  Array1dc nj, nk;
  Array1dc dnj, dnk;
  // Damping terms
  Array1dr nui, nue;
  // Thermal damping
  Array1dr nue_t;
  // Time spectral coefficients
  Array1dr ct, st;

  // Solver buffer variables with suffix
  // h   : variable one time step ahead
  // est : tilda variable one time step ahead
  Array1dc Ej2, Ek2, Eestk, Eestj, Ejh, Ekh, Ejh2, Ekh2;
  Array1dc nj2, nkh, njh;
  Array1dc nEj, nEk, nEestj, nEestk;

  Array1dc lDecay;
  Array1dc lDecay_t;
  Array1dr sDecay0, sDecay1, sDecay2, sDecay3, sDecay4;

  // Source terms variables
  Array1dc source_nk, source_Ek;
  Array1dc delta_nk, delta_Ek;

  // Variables for Box-Muller algorithm
  Array1dc iphi;
  Array1dr rho;
  //
  typedef ranlib::MersenneTwister MersenneTwister;
  typedef ranlib::sharedState sharedState; // shares an IRNG with other RNGs
  typedef ranlib::independentState independentState; // contains its own IRNG
  ranlib::UniformOpenClosed<Real, MersenneTwister, independentState> rngE;
  ranlib::UniformOpenClosed<Real, MersenneTwister, independentState> rngn;

  // Conserved quantities
  Real NN, PP, HH;

  // Pump value
  Complex Ep, Eph;

  // <n_jE_j> to compute SEE spectra
  Array1dc nEjmean;

  // Dissipated energy
  Real RR;

  // total electrostatic energy to thermal energy
  Real WW;

  // initial values of nue(k) and Fe(k) for diffusion
  Array1dr xd, nue0d, Fe0d;

  void initParsing (int nargs, char *args[]);
  void paramParsing();
  void checkParam() const;

  void zeroPadding(Array1dr &A) {
    A(k0) = Real(0.0);
  }
  void zeroPadding(Array1dc &A) {
    A(k0) = Complex(0.0, 0.0);
  }
  void initZeroPadding();
  void initCtSt();
  void allocateSolverBuffers();
  void allocateSourcesBuffers();

  void timeIntegrateClassic();
  void timeIntegratePump();
  void timeIntegrateFixedn();
  void timeIntegrateFixedE();
  void timeIntegrateLinear();
  void timeIntegrateNoFp();
  void (ZakharovSolver::*timeIntegrate[last_intg])();
  void calculateNPH();
  void calculateR();
  void calculateW();
  void calculateLangmuirDecay();
  void calculateThermalLangmuirDecay();
  void calculateSoundDecay();

  Array1dc getComplexGaussianNoise(const Array1dc & iphi, const Array1dr & rho,
                                   const Range & r=Range::all()) {
    return Array1dc(complexExp(iphi(r)) * sqrt(-2.0*log(rho(r))));
  }
  void calculateNkSource();
  void calculateEkSource();

  void printStatusInfo(std::ostream &os);

#if defined(HAVE_PLPLOT)
  void initPlot();
  void plot();
#endif

  struct simuTime {
    double ntime, utime;
    simuTime(Real _ntime, Real _utime) : ntime(_ntime), utime(_utime) {}
    friend std::ostream& operator<<(std::ostream &os, const simuTime &x) {
      std::ios::fmtflags f = os.flags() & std::ios::floatfield;
      int p = os.precision();
      return os << "Time    = "
             << std::fixed << std::setprecision(2) << std::setw(12)
             << timeutils::smartTime(x.utime)
             << " [" << x.ntime << " tau units]"
             << std::resetiosflags(std::ios::fixed)
             << std::setiosflags(f) << std::setprecision(p);
    }
  };
};

#endif // ZAKHAROV_H

