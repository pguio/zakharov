/**************************************************************************
 *
 * $Id: zakharov.cpp,v 1.233 2019/05/10 16:54:43 patrick Exp $
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

#include <zakharov.h>

#ifdef BZ_HAVE_COMPLEX_MATH
using blitz::conj;
#endif
using blitz::cos;
using blitz::exp;
using blitz::imag;
using blitz::log;
using blitz::max;
using blitz::pow2;
using blitz::real;
using blitz::sin;
using blitz::sqrt;
using blitz::where;
using blitz::zip;

using parser::header;
using parser::map_elt;
using parser::yesno;

using std::cout;
using std::endl;
using std::ios;
using std::ostream;

ostream& operator<<(ostream& os, const ZakharovSolver &s)
{

  os << header("1-d periodic Zakharov equation solver Version " VERSION, 1, 1)

     << "Initialisation mode       = " << map_elt(s.init_map,s.init_mode) << '\n'
     << "Time integration mode     = " << map_elt(s.intg_map,s.intg_mode) << '\n'

     << "##\n"
     << "L                         = " << s.L << '\n'
     << "h                         = " << s.h << '\n'
     // wave equation u_tt = u_xx
     << "CFL dt/dx                 = " << s.h/(s.L/s.N) << '\n'
     // wave equation u_t = i u_xx
     << "CFL dt/dx^2               = " << s.h/pow2(s.L/s.N) << '\n'
     << "FFT setup                 = " << s.fft << '\n'
#if defined(HAVE_FFTW3_FFT)
     << "Optimal plan for FFTW3    = " << yesno(s.optfftw3) << '\n'
#endif
     << "Number of iterations      = " << s.maxiter << '\n'

     << "##\n"
     << "Period of status display  = " << s.modulo_display << '\n'
#if defined(HAVE_PLPLOT)
     << "PLplot output             = " << yesno(s.plplot) << '\n'
#endif

     << "##\n"
     << "Density stochastic drive  = " << yesno(s.nk_source) << '\n'
     << "E-field stochastic drive  = " << yesno(s.Ek_source) << '\n'
     << "Stochastic drive seeds    = " << s.seed << '\n'

     << "##\n"
     << "Acoustic waves damping    = " << yesno(s.sound_damping) << '\n'
     << "Langmuir waves damping    = " << yesno(s.langmuir_damping) << '\n'

     << "##\n"
     << "Phase space diffusion     = " << yesno(s.diffusion) << '\n';

  if (s.langmuir_damping && blitz::min(s.nue) < 0)
    os << "##\n"
       << "Max Langmuir growth rate  = " << -blitz::min(s.nue)
       << " @ k = " << s.k(blitz::minIndex(s.nue)) << '\n';
  else
    os << "##\n"
       << "No Langmuir mode with growth rate\n";

  if (s.time_thermal > 0.0 || s.w_thermal > 0.0)
    os << "##\n";
  if (s.time_thermal > 0.0)
    os << "Time to thermal damping   = " << s.time_thermal << '\n';
  if (s.w_thermal > 0.0)
    os << "W to thermal damping      = " << s.w_thermal << '\n';

  os << s.mapper << '\n'
     << *s.mode << '\n'
     << s.spectra << '\n'
     << s.diagnostics;

  return os;
}

ZakharovSolver::ZakharovSolver(int nargs, char *args[]) : Parser(nargs, args),
  nk_source(false), Ek_source(false), seed(1),
  sound_damping(false), langmuir_damping(false),
  diffusion(false),
  time_thermal(-1), w_thermal(-1), thermal(false), optfftw3(false),
  iter(0), maxiter(0), modulo_display(100),
  mapper(nargs, args),
  init_mode(_solitons), intg_mode(_classic),
  spectra(nargs, args), diagnostics(nargs, args)
#if defined(HAVE_PLPLOT)
  ,plsnEx(0), plsnEk(0), plsdiffk(0), plplot(false)
#endif
{
  initParsing(nargs, args);

  Factory::instance().init();

  Factory::const_iterator i = Factory::instance().begin();
  Factory::const_iterator e = Factory::instance().end();
  for (; i!= e ; ++i) {
    const Factory::IDKeyType modeName(Factory::instance().getKey(i));
  }

  // initialisation mode map
  init_map.insert(Pair(_solitons, "solitons"));
  init_map.insert(Pair(_efluct  , "efluct"));
  init_map.insert(Pair(_tfluct  , "tfluct"));
  init_map.insert(Pair(_egauss  , "egauss"));
  init_map.insert(Pair(_densfluct  , "densfluct"));
  init_map.insert(Pair(_densparab  , "densparab"));

  // integration mode mode
  intg_map.insert(Pair(_classic, "classic"));
  intg_map.insert(Pair(_pump, "pump"));
  intg_map.insert(Pair(_fixed_n, "fixed density"));
  intg_map.insert(Pair(_fixed_E, "fixed E-field"));
  intg_map.insert(Pair(_linear, "linear"));
  intg_map.insert(Pair(_noFp, "no ponderomotive force"));

  timeIntegrate[_classic] = &ZakharovSolver::timeIntegrateClassic;
  timeIntegrate[_pump]    = &ZakharovSolver::timeIntegratePump;
  timeIntegrate[_fixed_n] = &ZakharovSolver::timeIntegrateFixedn;
  timeIntegrate[_fixed_E] = &ZakharovSolver::timeIntegrateFixedE;
  timeIntegrate[_linear]  = &ZakharovSolver::timeIntegrateLinear;
  timeIntegrate[_noFp]    = &ZakharovSolver::timeIntegrateNoFp;

  paramParsing();

  checkParam();

  mode = Factory::instance().create(init_map[init_mode], nargs, args);
}

ZakharovSolver::~ZakharovSolver()
{
  delete mode;
}

#define PARSE(FUNC)                 \
bool ZakharovSolver::FUNC() const   \
{                                   \
	if (Parser::FUNC()) {             \
		mapper.FUNC();                  \
		(*mode).FUNC();                 \
		spectra.FUNC();                 \
		diagnostics.FUNC();             \
		return true;                    \
	}                                 \
	return false;                     \
}                                   \
 

PARSE(parseHelp)

PARSE(parseVersion)

PARSE(parseTemplate)

void ZakharovSolver::initParsing(int nargs, char *args[])
{
  const string id(PACKAGE);
  const string version(
    VERSION "\n"
    "$Id: zakharov.cpp,v 1.233 2019/05/10 16:54:43 patrick Exp $"
    "\n");
  const string copyright(ZAK_COPYRIGHT);

  registerClass("ZakharovSolver");
  registerPackage(id, version, copyright);

  using parser::types::boolean;
  using parser::types::integer;
  using parser::types::real;
  using parser::types::intVect;

  insertOption(_init_mode  , "init_mode"   , integer, "Initialisation mode"       , Any(init_mode));

  insertOption(_intg_mode  , "intg_mode"   , integer, "Time integration mode"      , Any(intg_mode));

  insertOption(_maxiter, "maxiter", integer, "Number of time integration", Any(maxiter));
  insertOption(_modulo_display,"modulo_display", integer, "Period for status display", Any(modulo_display));

  insertOption(_nk_source, "nk_source", boolean, "Density fluctuation source", Any(nk_source));
  insertOption(_Ek_source, "Ek_source", boolean, "Electric field fluctuation source", Any(Ek_source));
  insertOption(_seed, "seed", intVect, "Seeds for stochastic drive RNGs (E,n)", Any(seed));

  insertOption(_sound_damping, "sound_damping", boolean, "Sound wave damping", Any(sound_damping));
  insertOption(_langmuir_damping, "langmuir_damping", boolean, "Langmuir wave damping", Any(langmuir_damping));

  insertOption(_diffusion, "diffusion", boolean, "Phase space diffusion", Any(diffusion));

  insertOption(_time_thermal, "time_thermal", real, "Time to switch to thermal [s]", Any(time_thermal));

  insertOption(_w_thermal, "w_thermal", real, "Wtotal to switch to thermal ", Any(time_thermal));

#if defined(HAVE_FFTW3_FFT)
  insertOption(_optfftw3, "optfftw3", boolean, "Optimal FFTW3 plan", Any(optfftw3));
#endif

#if defined(HAVE_PLPLOT)
  insertOption(_plplot, "plplot", boolean, "PLplot output", Any(plplot));
#endif
}

void ZakharovSolver::paramParsing ()
{

  parseOption(_init_mode, init_mode);
  parseOption(_intg_mode, intg_mode);

  parseOption(_maxiter, maxiter);
  parseOption(_modulo_display, modulo_display);

  parseOption(_nk_source, nk_source);
  parseOption(_Ek_source, Ek_source);
  if (nk_source || Ek_source ) {
    parseOption(_seed, seed);
    rngE.seed(seed(0));
    rngn.seed(seed(1));
  }

  parseOption(_sound_damping, sound_damping);
  parseOption(_langmuir_damping, langmuir_damping);

  parseOption(_diffusion, diffusion);

  parseOption(_time_thermal, time_thermal);
  parseOption(_w_thermal, w_thermal);

#if defined(HAVE_FFTW3_FFT)
  parseOption(_optfftw3, optfftw3);
#endif

#if defined(HAVE_PLPLOT)
  parseOption(_plplot, plplot);
#endif
}

void ZakharovSolver::checkParam() const
{
  checkMap(_init_mode, init_map, init_mode);
  checkMap(_intg_mode, intg_map, intg_mode);
}

void ZakharovSolver::initialise()
{
  mapper.initialise();

  N = mapper.getGridSize();
  M = mapper.getZeroPaddingSize();
  h = mapper.getNormalisedTimeStep();
  L = mapper.getNormalisedSystemLength();
  hOver2 = 0.5*h;
  ihOver2 = Complex(0.0,hOver2);

  fft.resize(N, -1);
#if defined(HAVE_FFTW3_FFT)
  if (optfftw3)
    fft.setPlanFlag(fourier::optimalPlanFlag);
#endif

#if defined(HAVE_PLPLOT)
  initPlot();
#endif

  initZeroPadding();
  xj.resize(N);
  xj = mapper.getNormalisedGridPoints();
  k.resize(N);
  k = mapper.getNormalisedWaveNumbers();

  k2.resize(N);
  k2 = pow2(k);

  ik2.resize(N);
  ik2 = zip(0.0, k2, Complex());

  hOver2k2.resize(N);
  hOver2k2 = hOver2 * k2;

  nui.resize(N);
  if (sound_damping) nui = mapper.getSoundWavesDamping();
  else nui = 0.0;

  nue.resize(N);
  nue_t.resize(N);
  if (langmuir_damping) {
    if (diffusion) {
      nue = mapper.initDiffusion();
      // keep copy of axis, nue and Fe at t=0
      xd.resize(2*M+1);
#if 0
      xd = getDiffusionAlphas();
#else
      Array1dr allkc(fourier::fftshift(k));
      xd = allkc(blitz::Range(N/2-M,N/2+M));
#endif
      nue0d.resize(2*M+1);
      nue0d = mapper.getDiffusionLangmuirWavesDamping();
      Fe0d.resize(2*M+1);
      Fe0d = mapper.getDiffusionElectronPdf();
    } else {
      nue = mapper.getLangmuirWavesDamping();
    }
    nue_t = mapper.getMaxwellianLangmuirWavesDamping();
  } else {
    nue = 0.0;
    nue_t = 0.0;
  }

  initCtSt();

  (*mode).initialise(mapper);
  (*mode).initEndn(mapper, fft, Ej, nj, dnj, Ek, nk, dnk);
  zeroPadding(Ek);
  fft.inverse(Ek , Ej);
  zeroPadding(nk);
  fft.inverse(nk , nj);
  zeroPadding(dnk);
  fft.inverse(dnk, dnj);

  if (nk_source || Ek_source) {
    source_nk.resize(N);
    source_nk = 0.0;

    source_Ek.resize(N);
    source_Ek = 0.0;

    (*mode).getSources(source_nk, source_Ek);
#if 0
    // Euler-Maruyuma
    source_nk *= sqrt(h)*sqrt(2.0*nui*pow2(k)/(pow2(k)+4.0*pow2(nui)));
#else
    // Exact Ornstein-Uhlenbeck process
    // See Gillespie D, Exact numerical simulation of the
    // Ornstein-Uhlenbeck process and its integral
    // Phys Rev E, 54, 2084-2091, 1996
    source_nk *= sqrt(1.0/(2.0*nui)*(1.0-exp(-2.0*nui*h)))*
                 sqrt(2.0*nui*pow2(k)/(pow2(k)+4.0*pow2(nui)));
#endif
    // set source to zero if no damping at k=0
    if (nui(0) == 0.0) source_nk(0) = 0.0;


#if 0
    // Euler-Maruyuma
    source_Ek *= sqrt(h)*sqrt(nue_t);
#else
    // Exact Ornstein-Uhlenbeck process
    source_Ek *= sqrt(1.0/(2.0*nue_t)*(1.0-exp(-2.0*nue_t*h)))*sqrt(nue_t);
#endif
    // set source to zero if no damping at k=0
    if (nue_t(0) == 0.0) source_Ek(0) = 0.0;

    zeroPadding(source_nk);
    zeroPadding(source_Ek);

    allocateSourcesBuffers();

    Range isaved(Range(0,adjfn(0,N-1,8),8));
    saveMatlab("sources.m", ios::out|ios::trunc, "k", k(isaved));
    saveMatlab("sources.m", ios::out|ios::app, "snk", source_nk(isaved));
    saveMatlab("sources.m", ios::out|ios::app, "sEk", source_Ek(isaved));
  }

  spectra.initialise(mapper);

  diagnostics.startInit();
  diagnostics.initnE(mapper);
  diagnostics.initNPH(mapper);
  NN = 0.0;
  PP = 0.0;
  HH = 0.0;
  diagnostics.initW(mapper);
  WW = 0.0;
  if (langmuir_damping) {
    diagnostics.initR(mapper);
    RR = 0.0;
  }
  if (diffusion) diagnostics.initDiffusion(mapper);
  diagnostics.initSk2(mapper, spectra);
  diagnostics.initSkf2(mapper, spectra);
  diagnostics.endInit();

  allocateSolverBuffers();

  calculateLangmuirDecay();
  calculateSoundDecay();
  calculateThermalLangmuirDecay();

  logfile = diagnostics.getDiagnosticFilename();
  logfile.append(".log");

  logs.open(logfile.c_str(), ios::out|ios::trunc);

  Parser::printCmd(logs);

  logs.setf(std::ios::scientific);
  logs << *this << endl;
}

void ZakharovSolver::calculateNPH()
{
  // from Payne et al., Numerical solution of the Zakharov equations,
  // J. Comp. Phys., 50, 482-498, 1983

  // Wrap around the last point
  Array1dr X(N+1);
  X(Range(0,N-1)) = xj;
  X(N) = -X(0);

  // Calculation of plasmon numbers (or action N) (Eq. 34)
  Array1dr Ni(N+1);
  Ni(Range(0,N-1)) = pow2(real(Ej))+pow2(imag(Ej));
  Ni(N) = Ni(0);
  NN = integrate1d(X, Ni);

#if 0
  {
    Real NNi, NNk;
    Array1dr Ek2(AbsSquare(Ek));
    Array1dr K(N+1), Nk(N+1);
    Nk = fourier::fftshift(Ek2);
    Nk(N) = Nk(0);
    K = blitz::tensor::i;
    // Parseval theorem for discrete Fourier transform (dimensionless)
    // sum |E_j|^2 = 1/N sum |E_k|^2
    NNi = integrate1d(K, Ni);
    NNk = 1./N*integrate1d(K, Nk);
    std::cout << "NN=" << NNi << " , " << NNk << std::endl;
    K(Range(0,N-1)) = fourier::fftshift(k);
    K(N) = -K(0);
    // Parseval theorem (in length/wavenumber physical units)
    // sum |E_j|^2 dx = sum |E_k/sqrt(2*pi)*L/N|^2 dk
    // where dx = L/N
    // and   dk = 2*pi/L
    NNi = integrate1d(X, Ni);
    NNk = pow2(L/N)/(2.0*m_pi)*integrate1d(K, Nk);
    std::cout << "NN=" << NNi << " , " << NNk << std::endl;
    // Parseval theorem (in time/frequency physical units)
    // sum |E_j|^2 dt = sum |E_k*tau|^2 df
    // where dt = tau
    // and   df = 1/(N*tau)
    // Parseval theorem (in time/angular frequency physical units)
    // sum |E_j|^2 dt = sum |E_k*tau/sqrt(2*pi)|^2 dw
    // where dt = tau
    // and   dw = 2*pi/(N*tau)
  }
#endif

  // calculation of the flux V (Eq. 37)
  Array1dc Vk(N);
  Vk = where(k==0.0, 0.0, -dnk/(I*k));
  Array1dc Vj(N);
  fft.inverse(Vk, Vj);

  // Calculation of the linear momentum P (Eq. 35)
  Array1dc ikEk(I*k*Ek);
  Array1dc dEj(N);
  fft.inverse(ikEk, dEj);

  Array1dc Pic(N);
#ifndef BZ_HAVE_COMPLEX_MATH
  Pic = I/2.0*(Ej*(real(dEj)-I*imag(dEj))-(real(Ej)-I*imag(Ej))*dEj)+nj*Vj;
#else
  Pic = I/2.0*(Ej*conj(dEj)-conj(Ej)*dEj)+nj*Vj;
#endif
  Array1dr Pi(N+1);
  Pi(Range(0,N-1)) = real(Pic);
  Pi(N) = Pi(0);
  PP = integrate1d(X, Pi);

  // Calculation of the energy H (Eq. 36)
  Array1dr Hi(N+1);
  Hi(Range(0,N-1)) = pow2(real(dEj))+pow2(imag(dEj))+
                     real(nj)*(pow2(real(Ej))+pow2(imag(Ej)))+
                     0.5*(pow2(real(nj))+pow2(real(Vj)));
  Hi(N) = Hi(0);
  HH = integrate1d(X, Hi);
}

void ZakharovSolver::calculateR()
{
  // R is defined in Hanssen et al., Numerical test of the weak turbulence
  // approximation to ionospheric Langmuir turbulence, JGR, 97, 12073-12091,
  // 1992. (Eqs. 15-16 p 12078)
  //
  Real nuec = mapper.getNormalisedElectronCollisionFreq();

  // total amount of Langmuir energy dissipated
  Real DL = 2.0*blitz::sum(where(nue >= 0.0,
                                 nue*(pow2(real(Ek))+pow2(imag(Ek))), 0.0));

  // energy dissipated by collisions
  Real DC = nuec*blitz::sum(where(nue >= 0.0,
                                  pow2(real(Ek))+pow2(imag(Ek)), 0.0));

  // ratio of energy dissipated by collisions to the total amount of
  // dissipation.
  // R = 1 all the dissipitation is collisional
  // R = 0 dissipation is completely non-collisional (Landau damping)
  RR = DC / DL;
}

void ZakharovSolver::calculateW()
{
  Real N2W = mapper.getN2Wconstant();
  WW = NN*N2W;

#if 0
  Real sk = sum(pow2(real(Ek))+pow2(imag(Ek)));
  Real sj = sum(pow2(real(Ej))+pow2(imag(Ej)));
  Real Erms = sum(pow2(real(Ej))+pow2(imag(Ej)))/N;
  std::cout
      << "1/N sum |E_k|^2=" << sk/N << ' '
      << "sum |E_j|^2=" << sj << '\n';
  std::cout << "N=sj*dx=" << sj*(xj(1)-xj(0)) << std::endl;
  std::cout << "N2W=" << N2W << std::endl;
  std::cout << "W=" << WW << " " << Erms*L*N2W << std::endl;
#endif
}

void ZakharovSolver::solve()
{
  if (! maxiter ) return;

  logs.setf(std::ios::scientific);

  timer.start();

  do {

    switchToThermal(logs);

    if (! (iter % modulo_display) || iter == maxiter) {
      printStatusInfo(cout);
      printStatusInfo(logs);
      cout << "nEjmean = " << nEjmean(0) << " " << blitz::sum(nj*Ej)/Real(N) << endl;
      logs << "nEjmean = " << nEjmean(0) << " " << blitz::sum(nj*Ej)/Real(N) << endl;
    }

#if defined(HAVE_PLPLOT)
    plot();
#endif

    if (diagnostics.saveStandardDiagnostics(iter)) {
      calculateNPH();
      diagnostics.savenE(mapper, iter, nj, Ej, nk, Ek);
      diagnostics.saveNPH(mapper, iter, NN, PP, HH);
      calculateW();
      diagnostics.saveW(mapper, iter, WW);
      if (langmuir_damping) {
        calculateR();
        diagnostics.saveR(mapper, iter, RR);
      }
      if (diffusion) diagnostics.saveDiffusion(mapper, iter);
    }

    spectra.updateSpectra(mapper, iter, nk, Ek, nEjmean);
    diagnostics.saveSk2(mapper, spectra);
    diagnostics.saveSkf2(mapper, iter, false, spectra);

    (this->*timeIntegrate[intg_mode])();

  } while (iter <= maxiter);

  spectra.updateSpectra(mapper, iter, nk, Ek, nEjmean);
  diagnostics.saveSk2(mapper, spectra);
  diagnostics.saveSkf2(mapper, iter, true, spectra);

  cout << endl;
}

void ZakharovSolver::switchToThermal(std::ostream &os)
{
  if (langmuir_damping && !thermal && (
        (time_thermal > 0 && time_thermal <= mapper.toPhysicalTime(h*iter)) ||
        (w_thermal > 0 && w_thermal <= WW)
      ) ) {

    cout  << "**********************************\n"
          << "Switched to thermal landau damping\n"
          << simuTime(h*iter, mapper.toPhysicalTime(h*iter)) << '\n'
          << "W       = " << WW << '\n'
          << "**********************************" << endl;

    logs  << "**********************************\n"
          << "Switched to thermal landau damping\n"
          << simuTime(h*iter, mapper.toPhysicalTime(h*iter)) << '\n'
          << "W       = " << WW << '\n'
          << "**********************************" << endl;

    lDecay = lDecay_t;
    thermal = true;
  }
}

void ZakharovSolver::initZeroPadding()
{
  using blitz::tensor::i;

  // symmetric or zero-centered indices
  // [-N/2:1:N/2-1] (even N),
  // [-N/2:1:N/2] (odd N)
  Array1di signedIndex(N);
  signedIndex = i - (N>>1);
  // shifted or FFT indices
  // [0:N/2-1 -N/2:-1:-1] (even N)
  // [0:N/2 -N/2:-1:-1] (odd N)
  signedIndex = fourier::ifftshift(signedIndex);
  // zero padding indices when
  // signedIndex > M or signedIndex < M
  // in shifted indices
  k0.setRange(M+1,(N-1)-M); // length N-1-2M k=M+1..N-1-M

  // indices of wavenumbers k>=0 (kp) and k<0 (kn) to be computed,
  // i.e. within the zero padding limits
  kp.setRange(0, M);    // length M+1 : k=0 k=1..M
  kn.setRange(N-M,N-1); // length M   : k=N-1..N-M
}

void ZakharovSolver::initCtSt()
{
  using blitz::sinh;
  using blitz::cosh;

  // Calculation of the cos and sin terms in Eq. 10
  Array1dr k2_nui2(k2 - pow2(nui));

  // cos(i x) = cosh(x)
  ct.resize(N);
  ct = where(k2_nui2 > 0.0, cos(sqrt(k2_nui2)*h), cosh(sqrt(-(k2_nui2))*h));

  // sin(i x) = i sinh(x)
  st.resize(N);
  st = where(k2_nui2 > 0.0,
             sin(sqrt(k2_nui2)*h)/sqrt(k2_nui2), where(k2_nui2 != 0.0,
                 sinh(sqrt(-(k2_nui2))*h)/sqrt(-(k2_nui2)), h));
}

void ZakharovSolver::allocateSolverBuffers()
{
  Ej2.resize(N);
  Ek2.resize(N);
  Ek2 = Complex(0.0,0.0);

  Eestj.resize(N);
  Eestk.resize(N);
  Eestk = Complex(0.0,0.0);

  Ejh.resize(N);
  Ekh.resize(N);
  Ekh = Complex(0.0,0.0);

  Ejh2.resize(N);
  Ekh2.resize(N);
  Ekh2 = Complex(0.0,0.0);

  nj2.resize(N);

  njh.resize(N);
  nkh.resize(N);
  nkh = Complex(0.0,0.0);

  nEj.resize(N);
  nEk.resize(N);
  nEk = Complex(0.0,0.0);
  nEestj.resize(N);
  nEestk.resize(N);
  nEestk = Complex(0.0,0.0);

  lDecay.resize(N);
  lDecay_t.resize(N);
  sDecay0.resize(N);
  sDecay1.resize(N);
  sDecay2.resize(N);
  sDecay3.resize(N);
  sDecay4.resize(N);

  // Array to contain <n_jE_j> to compute SEE frequency spectra
  nEjmean.resize(1);
  nEjmean = Complex(0.0,0.0);
}

void ZakharovSolver::allocateSourcesBuffers()
{
  delta_nk.resize(N);
  delta_Ek.resize(N);

  iphi.resize(N);
  iphi = Complex(0.0,0.0);

  rho.resize(N);
  rho = 1.0;
}

void ZakharovSolver::calculateLangmuirDecay()
{
  // Calculate the term exp(-(ik2+nue)*h) in Eq. 11
  Array1dc z(-(ik2+nue)*h);
  lDecay = complexExp(Array1dc(-(ik2+nue)*h));

  Range isaved(Range(0,adjfn(0,N-1,8),8));
  saveMatlab("ldecay.m", ios::out|ios::trunc, "k", k(isaved));
  saveMatlab("ldecay.m", ios::out|ios::app, "ik2", ik2(isaved));
  saveMatlab("ldecay.m", ios::out|ios::app, "ld", lDecay(isaved));
}

void ZakharovSolver::calculateThermalLangmuirDecay()
{
  // Calculate the term exp(-(ik2+nue)*h) in Eq. 11
  Array1dc z(-(ik2+nue_t)*h);
  lDecay_t = complexExp(Array1dc(-(ik2+nue_t)*h));

  Range isaved(Range(0,adjfn(0,N-1,8),8));
  saveMatlab("ldecayt.m", ios::out|ios::trunc, "k", k(isaved));
  saveMatlab("ldecayt.m", ios::out|ios::app, "ik2", ik2(isaved));
  saveMatlab("ldecayt.m", ios::out|ios::app, "ld", lDecay_t(isaved));
}

void ZakharovSolver::calculateSoundDecay()
{

  sDecay0 = exp(-nui*h);
  sDecay1 = sDecay0*(ct+nui*st); // nk  term in Eq. 10
  sDecay2 = sDecay0*st;          // dnk term in Eq. 10
  sDecay3 = hOver2k2*sDecay0*st; // Ek2 term in Eq. 10
  sDecay4 = exp(-2.0*nui*h);

  Range isaved(Range(0,adjfn(0,N-1,8),8));
  saveMatlab("sdecay.m", ios::out|ios::trunc, "k", k(isaved));
  saveMatlab("sdecay.m", ios::out|ios::app, "ct", ct(isaved));
  saveMatlab("sdecay.m", ios::out|ios::app, "st", st(isaved));
  saveMatlab("sdecay.m", ios::out|ios::app, "sd0", sDecay0(isaved));
  saveMatlab("sdecay.m", ios::out|ios::app, "sd1", sDecay1(isaved));
  saveMatlab("sdecay.m", ios::out|ios::app, "sd2", sDecay2(isaved));
  saveMatlab("sdecay.m", ios::out|ios::app, "sd3", sDecay3(isaved));
  saveMatlab("sdecay.m", ios::out|ios::app, "sd4", sDecay4(isaved));

}


void ZakharovSolver::calculateNkSource()
{
  if (nk_source) {

    iphi(0) = Complex(0.0, m_2pi * rngn.random());
    rho(0)  = rngn.random();
    for (int i=1; i<=M; ++i) {
      // Make the density spectral fluctuation such that
      // nk( -k ) = conj( nk( k ) )
      iphi(i)   = Complex(0.0, m_2pi * rngn.random());
      iphi(N-i) = -iphi(i);
      rho(i)    = rngn.random();
      rho(N-i)  = rho(i);
    }
    delta_nk(kp) = getComplexGaussianNoise(iphi, rho, kp)*source_nk(kp);
    delta_nk(kn) = getComplexGaussianNoise(iphi, rho, kn)*source_nk(kn);
    delta_nk(0) = Complex(real(delta_nk(0)), 0.0);

    nkh += delta_nk;

    if (! iter) {
      Range isaved(Range(0,adjfn(0,N-1,8),8));
      saveMatlab("delta_nk.m", ios::out|ios::trunc, "k", k(isaved));
      saveMatlab("delta_nk.m", ios::out|ios::app, "dnk", delta_nk(isaved));
    }
  }
}

void ZakharovSolver::calculateEkSource()
{
  if (Ek_source) {

    iphi(0) = Complex(0.0, m_2pi * rngE.random());
    rho(0)  = rngE.random();
    for (int i=1; i<=M; ++i) {
      iphi(i)   = Complex(0.0, m_2pi * rngE.random());
      iphi(N-i) = Complex(0.0, m_2pi * rngE.random());
      rho(i) = rngE.random();
      rho(N-i) = rngE.random();
    }

    delta_Ek(kp) = getComplexGaussianNoise(iphi, rho, kp)*source_Ek(kp);
    delta_Ek(kn) = getComplexGaussianNoise(iphi, rho, kn)*source_Ek(kn);

    Ekh += delta_Ek;

    if (! iter) {
      Range isaved(Range(0,adjfn(0,N-1,8),8));
      saveMatlab("delta_Ek.m", ios::out|ios::trunc, "k", k(isaved));
      saveMatlab("delta_Ek.m", ios::out|ios::app, "dEk", delta_Ek(isaved));
    }
  }
}

void ZakharovSolver::printStatusInfo(std::ostream &os)
{
  os << "Iteration " << timeutils::iterStatus(iter,maxiter);
  if (iter > 0) {
    double usedTime = timer.realRunningElapsed();
    double remTime = usedTime*(maxiter-iter)/iter;
    os
        << " (used time " << timeutils::toHMS(usedTime)
        << ", remaining " << timeutils::toHMS(remTime) << ")" << endl;
  } else {
    os << '\n';
  }

  os << simuTime(h*iter, mapper.toPhysicalTime(h*iter)) << endl;

  calculateNPH();

  os << "N, P, H = " << NN << ", " << PP << ", " << HH << endl;

  if (langmuir_damping) {
    calculateR();
    os << "R       = " << RR << endl;
  }

  calculateW();
  os << "W       = " << WW << endl;

}


#if defined(HAVE_PLPLOT)

void ZakharovSolver::initPlot()
{
  if (plplot) {
    // general plot n(x), E(x)
    plsnEx = new plstream();
    plsnEx->sdev("xwin");
    plsnEx->init();           // Initialize PLplot
    plsnEx->fontld(1);        // Select the multi-stroke font
    plsnEx->font(1);          // Select font 1
    plsnEx->spause(false);    // Turn off pause to make this a slave
    // general plot n(k), E(k)
    plsnEk = new plstream();
    plsnEk->sdev("xwin");
    plsnEk->init();           // Initialize PLplot
    plsnEk->fontld(1);        // Select the multi-stroke font
    plsnEk->font(1);          //
    plsnEk->spause(false);
    // nue and Fe
    if (diffusion) {
      plsdiffk = new plstream();
      plsdiffk->sdev("xwin");
      plsdiffk->init();           // Initialize PLplot
      plsdiffk->fontld(1);        // Select the multi-stroke font
      plsdiffk->font(1);          //
      plsdiffk->spause(false);
    }
    // max size is FFT size
    xp = new PLFLT[N];
    yp = new PLFLT[N];
  }
}


void ZakharovSolver::plot()
{
  if (plsnEx && (! (iter % modulo_display) || iter == maxiter) ) {
    plsnEx->eop();
    plsnEx->bop();
    plsnEx->ssub(1,2);
    for (int i=0; i < N; ++i) {
      xp[i] = xj(i);
      yp[i] = real(Ej(i));
    }
    xpmin = min(xj);
    xpmax = max(xj);
    ypmin = -max(abs(Ej));
    ypmax = max(abs(Ej));
    plsnEx->col0(15); // white
    plsnEx->env( xpmin, xpmax, ypmin, ypmax, 0, 0 );
    plsnEx->col0(14); // salmon
    plsnEx->line(N, xp, yp);

    for (int i=0; i < N; ++i) {
      yp[i] = imag(Ej(i));
    }
    plsnEx->col0(13); // magenta
    plsnEx->line(N, xp, yp);

    for (int i=0; i < N; ++i) {
      yp[i] = abs(Ej(i));
    }
    plsnEx->col0(12); // turquoise
    plsnEx->line(N, xp, yp);

    plsnEx->col0(15);
    ostringstream os;
    os << "t = " << h*iter;
    plsnEx->lab( "x", "E(x)", os.str().c_str() );

    for (int i=0; i < N; ++i) {
      yp[i] = real(nj(i));
    }
    ypmin = min(real(nj));
    ypmax = max(real(nj));
    plsnEx->col0(15);
    plsnEx->env( xpmin, xpmax, ypmin, ypmax, 0, 0 );
    plsnEx->col0(14);
    plsnEx->line(N, xp, yp);

    plsnEx->col0(15);
    plsnEx->lab( "x", "n(x)", "");
  }

  if (plsnEk && (! (iter % modulo_display) || iter == maxiter) ) {
    plsnEk->eop();
    plsnEk->bop();
    plsnEk->ssub(1,2);
    int Np = 2*M+1;
    Array1dr tmpr(fourier::fftshift(k));
    Array1dr kc(Np);
    kc = tmpr(blitz::Range(N/2-M,N/2+M));
    Array1dc tmpc(fourier::fftshift(Ek));
    Array1dc Ekc(Np);
    Ekc = tmpc(blitz::Range(N/2-M,N/2+M));
    Real mx = max(abs(Ekc));
    for (int i=0; i < Np; ++i) {
      xp[i] = kc(i);
      //yp[i] = abs(Ekc(i);
      yp[i] = 20.0 * log10(abs(Ekc(i))/mx);
    }
    xpmin = min(kc);
    xpmax = max(kc);
    //ypmin = 0;
    //ypmax = max(abs(Ekc));
    ypmin = min(20.0 * log10(abs(Ekc)/mx));
    ypmax = 0;
    plsnEk->col0(15);
    plsnEk->env( xpmin, xpmax, ypmin, ypmax, 0, 0 );
    plsnEk->col0(14);
    plsnEk->line(Np, xp, yp);

    plsnEk->col0(15);
    ostringstream os;
    os << "t = " << h*iter;
    plsnEk->lab( "k", "|E(k)|", os.str().c_str() );

    tmpc = fourier::fftshift(nk);
    Array1dc nkc(Np);
    nkc = tmpc(blitz::Range(N/2-M,N/2+M));
    mx = max(abs(nkc));
    for (int i=0; i < Np; ++i) {
      //yp[i] = abs(nkc(i));
      yp[i] = 20.0 * log10(abs(nkc(i))/mx);
    }
    //ypmin = 0;
    //ypmax = max(abs(nkc));
    ypmin = min(20.0 * log10(abs(where(kc == 0.0, 1.0, nkc))/mx));
    ypmax = 0;
    plsnEk->col0(15);
    plsnEk->env( xpmin, xpmax, ypmin, ypmax, 0, 0 );
    plsnEk->col0(14);
    plsnEk->line(Np, xp, yp);

    plsnEk->col0(15);
    plsnEk->lab( "k", "|n(k)|", "");

    plsnEk->flush();
  }

  if (diffusion && plsdiffk && (! (iter % modulo_display) || iter == maxiter) ) {
    plsdiffk->eop();
    plsdiffk->bop();
    plsdiffk->ssub(1,2);
    int Np = 2*M+1;
    Array1dr nuediff(mapper.getDiffusionLangmuirWavesDamping());
    Array1dr Fediff(mapper.getDiffusionElectronPdf());
    for (int i=0; i < Np; ++i) {
      xp[i] = xd(i);
      yp[i] = nuediff(i);
    }
    blitz::TinyVector<int,1> imin = blitz::minIndex(nue0d); // max growth rate
#if 0
    // trick to get sign
    // http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
    //int sign = (xd(imin) > 0) - (xd(imin) < 0);
    xpmin = xd(imin) - 4.0 * abs(xd(imin)); // min(xd);
    xpmax = xd(imin) + 4.0 * abs(xd(imin)); // max(xd);
#else
    xpmin = 0.5 * min(xd);
    xpmax = 0.5 * max(xd);
#endif
    ypmin = 2.0 * min(nue0d);
    ypmax = 3.0 * abs(min(nue0d));
    plsdiffk->col0(15);
    plsdiffk->env( xpmin, xpmax, ypmin, ypmax, 0, 0 );
    plsdiffk->col0(14);
    plsdiffk->line(Np, xp, yp);

    plsdiffk->col0(15);
    ostringstream os;
    os << "t = " << h*iter;
    plsdiffk->lab( "k", "#gn#de#u(k)", os.str().c_str() );

    for (int i=0; i < Np; ++i) {
      yp[i] = Fediff(i);
    }
    ypmin = 0;
    ypmax = 2.0 * Fediff(imin);
    plsdiffk->col0(15);
    plsdiffk->env( xpmin, xpmax, ypmin, ypmax, 0, 0 );
    plsdiffk->col0(14);
    plsdiffk->line(Np, xp, yp);

    plsdiffk->col0(15);
    plsdiffk->lab( "k", "F#de#u(k)", "");

    plsdiffk->flush();
  }
}
#endif
