/**************************************************************************
 *
 * $Id: mapper.cpp,v 1.121 2014/10/21 10:25:44 patrick Exp $
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

#include <mapper.h>

#if defined(HAVE_GNUPLOT)
#include <gnuplot-interface.h>
#endif


using blitz::abs;
using blitz::exp;
using blitz::max;
using blitz::pow2;
using blitz::pow3;
using blitz::pow4;
using blitz::sqrt;
using blitz::where;

using blitz::imag;
using blitz::real;

using parser::header;

using std::imag;
using std::real;

using std::cout;
using std::endl;
using std::ios;
using std::ostream;

ostream& operator<<(ostream& os, const Mapper &m)
{
  static const Real Mega=1e6;
  os
      << header("Mapper setup")

      << "N        = " << m.N << '\n'
      << "M        = " << m.M << '\n';

  if (m.physical_units) {
    os
        << "L        = " << m.L_m << " [m]\n"
        << "h        = " << m.h_s << " [s]\n"

        << "##\n"
        << "ne       = " << m.ne << " [m-3]\n"
        << "Te       = " << m.Te << " [K]\n"
        << "nuec     = " << m.nuec << " [s-1]\n"
        << "Ne       = " << m.Ne << '\n'

        << "##\n"
        << "mi       = " << m.mni << " [amu]\n"
        << "Ti       = " << m.Ti << " [K]\n"
        << "nuic     = " << m.nuic << " [s-1]\n";

    if (!m.nb_ne.empty()) {
      os
          << "##\n"
          << "nb/ne    = " << m.nb_ne << '\n'
          << "ub       = " << m.ub << " [m s-1]\n"
          << "dub/ub   = " << m.dub_ub << '\n'
          << "Eb       = " << m.Eb << " [eV]\n"
          << "kb       = " << m.kb << " [m-1]\n"
          << "Gamma_b  = " << m.Gam << '\n'
          << "W_L      = " << m.W_L << '\n'
          << "W_L cold = " << m.W_Lc << '\n';
    }

    if (m.Epump != 0.0) {
      os
          << "##\n"
          << "wpump    = " << m.wpump << " [rad s-1]\n"
          << "fpump    = " << m.wpump/m_2pi/Mega << " [MHz]\n"
          << "dOmega   = " << m.dOmega << " [rad s-1]\n"
          << "dOmega   = " << m.dOmega_n << "\n"
          << "Epump    = " << m.Epump << " [V m-1]\n"
          << "Epump    = " << m.Epump_n << "\n"
          << "kpump    = " << m.kpump << " [m-1]\n"
          << "kpump    = " << m.kpump_n << "\n";
    }

    if (m.vphi.size() > 0) {
      os
          << "##\n"
          << "tau_1/2  = " << 1.0/m.nb_ne[0]/m.wpe << " [s]\n"
          << "tau_1/2  = " << 1.0/m.nb_ne[0]*2.0/3.0*m.eta*SI::me/m.mi << '\n';
    }

    os
        << "##\n"
        << "eta      = " << m.eta << '\n'
        << "wpe      = " << m.wpe << " [rad s-1]\n"
        << "fpe      = " << m.wpe/m_2pi/Mega << " [MHz]\n"
        << "ve       = " << m.ve << " [m s-1]\n"
        << "vi       = " << m.vi << " [m s-1]\n"
        << "cs       = " << m.cs << " [m s-1]\n"
        << "lambdae  = " << m.lambdae << " [m]\n"

        << "##\n"
        << "tau      = " << m.tau << " [s]\n"
        << "chi      = " << m.chi << " [m]\n"
        << "epsilon  = " << m.epsilon << " [V m-1]\n"
        << "nu       = " << m.nu << " [m-3]";
  } else {
    os
        << "L        = " << m.L << '\n'
        << "h        = " << m.h;
  }
  return os;
}


Mapper::Mapper(int nargs, char *args[]) : Parser(nargs, args),
  N(NFFT_DEFAULT), M(NFFT_DEFAULT/3),
  physical_units(false),
  ne(1.0e11), Te(2000.0), nuec(0.0),
  nb_ne(0), ub(0), dub_ub(0), Eb(0), kb(0), Gam(0), W_L(0), W_Lc(0),
  mni(16), Ti(1000.0), nuic(0.0), dOmega(0.0), Epump(0.0),
  tau(1.0), chi(1.0), epsilon(1.0), nu(1.0)
#if 0 // buffer to estimate time mean value of <E(k)^2>
  , nt(20), it(0), bufferComplete(false)
#endif
{
  initParsing(nargs, args);
  paramParsing();
  checkParam();
}

Mapper::~Mapper()
{}

void Mapper::initParsing(int nargs, char *args[])
{
  const string id(PACKAGE);
  const string version(
    VERSION "\n"
    "$Id: mapper.cpp,v 1.121 2014/10/21 10:25:44 patrick Exp $"
    "\n");
  const string copyright(ZAK_COPYRIGHT);

  registerClass("Mapper");
  registerPackage(id, version, copyright);

  using parser::types::boolean;
  using parser::types::integer;
  using parser::types::real;
  using parser::types::realVect;


  insertOption(_physical_units, "physical_units", boolean, "Parameters in physical units", Any(physical_units));
  parseOption(_physical_units, physical_units);

  insertOption(__N, "N", integer, "Number of grid points"  , Any(N));
  insertOption(__M, "M", integer, "k-mode truncation index", Any(M));

  if (physical_units) {
    insertOption(__L, "L", real, "Spatial period         [m]", Any(L_m));
    insertOption(__h, "h", real, "Time increment         [s]", Any(h_s));

    insertOption(_ne, "ne", real, "Electron density     [m-3]", Any(ne));
    insertOption(_Te, "Te", real, "Electron temperature [K]", Any(Te));
    insertOption(_nuec, "nuec", real, "Electron collisions freq [s-1]", Any(nuec));
    insertOption(_Ne, "Ne", real, "Number of electrons", Any(Ne));

    insertOption(_nb_ne, "nb_ne", realVect, "Beam/background density ratio", Any());
    insertOption(_ub, "ub", realVect, "Beam velocity [m s-1]", Any());
    insertOption(_dub_ub, "dub_ub", realVect, "Ratio dub/ub", Any());

    insertOption(_mi, "mi", real, "Ion mass             [amu]", Any(mni));
    insertOption(_Ti, "Ti", real, "Ion temperature      [K]", Any(Ti));
    insertOption(_nuic, "nuic", real, "Ion collisions freq [s-1]", Any(nuic));

    insertOption(_dOmega, "dOmega", real, "Pump frequency mismatch [rad-1]", Any(dOmega));
    insertOption(_Epump, "Epump", real, "Pump amplitude [V m-1]", Any(Epump));

  } else {
    insertOption(__L, "L", real, "Normalised spatial period", Any(L));
    insertOption(__h, "h", real, "Normalised time increment", Any(h));
  }
}

void Mapper::paramParsing()
{
  parseOption(__N, N);
  M = N/3;
  parseOption(__M, M);

  if (physical_units) {
    parseOption(__L, L_m);
    parseOption(__h, h_s);

    parseOption(_ne, ne);
    parseOption(_Te, Te);
    parseOption(_nuec, nuec);
    Ne = ne*pow3(L_m);
    parseOption(_Ne, Ne);

    parseOption(_nb_ne, nb_ne);
    parseOption(_ub, ub);
    parseOption(_dub_ub, dub_ub);

    parseOption(_mi, mni);
    parseOption(_Ti, Ti);
    parseOption(_nuic, nuic);

    parseOption(_dOmega, dOmega);
    parseOption(_Epump, Epump);

  } else {
    parseOption(__L, L);
    parseOption(__h, h);
  }
}

void Mapper::checkParam() const
{
  if (nb_ne.size() != ub.size() || nb_ne.size() != dub_ub.size()) {
    throw ClassException("Mapper",
                         "nb_ne, ub and dub_ub variables should have same size");
  }
}

void Mapper::initialise()
{
  if (physical_units) {

    // mass number into mass
    mi = mni * SI::amu;

    // Specific heats ratio for electron/ions
    // 1 for isothermal electrons
    // 3 for adiabatic responding ions
    Real ce = 1.0;
    Real ci = 3.0;
    eta = (ce*Te+ci*Ti)/Te;

    // Plasma parameters
    wpe = sqrt(ne*pow2(SI::e)/(SI::eps0*SI::me));
    ve = sqrt(SI::Kb*Te/SI::me);
    vi = sqrt(SI::Kb*Ti/mi);
    cs = sqrt(SI::Kb*eta*Te/mi);
    ve2 = SI::Kb*Te/SI::me;
    lambdae = sqrt(SI::eps0*SI::Kb*Te/(ne*pow2(SI::e)));
    ke = m_2pi / lambdae;

    // Normalisation constants from
    // Nicholson, Introduction to plasma theory, 1992, p. 180
    // Forme, Parametric decay of beam-driven Langmuir wave and
    // enhanced ion-acoustic fluctuations, Ann. Geo., 17, 1172-1181, 1999
    // doc/zakharov_sources.tex
    tau = (3.0/2.0)*(mi/(SI::me*eta))*(1.0/wpe);
    chi = (3.0/2.0)*sqrt(mi/(SI::me*eta))*lambdae;
    epsilon = 4.0/sqrt(3.0)*eta*sqrt(SI::me/mi)*sqrt(ne*SI::Kb*Te/SI::eps0);
    nu = (4.0/3.0)*(eta*SI::me)/mi*ne;

    L = toNormalisedLength(L_m);
    h = toNormalisedTime(h_s);

    // Pump parameters
    wpump = dOmega + wpe;

    dOmega_n = toNormalisedFreq(dOmega);
    Epump_n = toNormalisedElectricField(Epump);
    // Eq. 7 in Hanssen et al., 1992
    kpump_n = 0.5*(sqrt(4.0*dOmega_n+1)-1.0);
    kpump = toPhysicalWaveNumber(kpump_n);

    // beam parameters
    if (nb_ne.size() > 0) {
      Eb.resize(nb_ne.size());
      kb.resize(nb_ne.size());
      Gam.resize(nb_ne.size());
      W_L.resize(nb_ne.size());
      W_Lc.resize(nb_ne.size());
      for (unsigned i=0; i<nb_ne.size(); ++i) {
        Eb[i] = (1.0/2.0)*SI::me*pow2(ub[i])/SI::eV;
        // wavenumber of the excited Langmuir
        kb[i] = wpe/ub[i];
        // Parameter as defined in Eq.14 in Guio and Forme, PoP, 2006
        Gam[i] = nb_ne[i]/dub_ub[i];
        // Time asymptotic energy in the electron beam instability
        W_L[i] = SI::me*nb_ne[i]*ne/(ub[i]-dub_ub[i]*ub[i])*
                 (
                   (1.0/3.0)*(pow3(ub[i])-pow3(dub_ub[i]*ub[i]))-
                   (1.0/2.0)*dub_ub[i]*ub[i]*(pow2(ub[i])-pow2(dub_ub[i]*ub[i]))
                 )/(ne*SI::Kb*Te);
        // Time asymptotic energy in the cold electron beam instability
        // 2/3 beam energy 1/2*me*nb*vb^2
        W_Lc[i] = (1.0/3.0)*SI::me*nb_ne[i]*ne*pow2(ub[i])/(ne*SI::Kb*Te);
      }
    }
    // primary background plasma parameters
    saveMatlab("phys.m", ios::out|ios::trunc, "ne", ne);
    saveMatlab("phys.m", ios::out|ios::app, "Te", Te);
    saveMatlab("phys.m", ios::out|ios::app, "nuec", nuec);
    saveMatlab("phys.m", ios::out|ios::app, "mi", mni);
    saveMatlab("phys.m", ios::out|ios::app, "Ti", Ti);
    saveMatlab("phys.m", ios::out|ios::app, "nuic", nuic);

    // derived background plasma parameters
    saveMatlab("phys.m", ios::out|ios::app, "eta", eta);
    saveMatlab("phys.m", ios::out|ios::app, "wpe", wpe);
    saveMatlab("phys.m", ios::out|ios::app, "ve", ve);
    saveMatlab("phys.m", ios::out|ios::app, "vi", vi);
    saveMatlab("phys.m", ios::out|ios::app, "cs", cs);
    saveMatlab("phys.m", ios::out|ios::app, "lambdae", lambdae);

    // scaling factors
    saveMatlab("phys.m", ios::out|ios::app, "tau", tau);
    saveMatlab("phys.m", ios::out|ios::app, "chi", chi);
    saveMatlab("phys.m", ios::out|ios::app, "epsilon", epsilon);
    saveMatlab("phys.m", ios::out|ios::app, "nu", nu);
  }
}

Array1dr Mapper::getPhysicalWaveNumbers() const
{
  using blitz::tensor::i;
  Array1dr k(N);
  // zero-centred
  k = m_2pi * (i - (N>>1))/toPhysicalLength(L);
  return k;
}

Array1dr Mapper::getNormalisedWaveNumbers() const
{
  using blitz::tensor::i;
  Array1dr k(N);
  // zero-centred
  k = m_2pi * (i - (N>>1))/L;
  // shifted with zero-component first
  k = fourier::ifftshift(k);
  return k;
}

Array1dr Mapper::getAlphas() const
{
  // zero-centred alpha defined as k \lambda_e
  return Array1dr( getPhysicalWaveNumbers() * lambdae);
}

Array1dr Mapper::getPhysicalGridPoints() const
{
  using blitz::tensor::i;
  Array1dr x(N);
  x = toPhysicalLength(L)*(1.0/N*i - 0.5);
  return x;
}

Array1dr Mapper::getNormalisedGridPoints() const
{
  using blitz::tensor::i;
  Array1dr x(N);
  x = L*(1.0/N*i - 0.5);
  return x;
}

Array1dr Mapper::getDiffusionAlphas() const
{
  Array1dr alphas(getAlphas());
  return Array1dr( alphas(blitz::Range(N/2-M,N/2+M)) );
}


Array1dr Mapper::getDiffusionElectronPdf() const
{
  return Array1dr( Fe(blitz::Range(N/2-M,N/2+M)) );
}

Array1dr Mapper::getDiffusionLangmuirWavesDamping() const
{
  return Array1dr( 0.5*nuec +
                   toPhysicalFreq(gamma(blitz::Range(N/2-M,N/2+M))) );
}

Array1dr Mapper::getnkFluctuationsLevel() const
{
  Array1dr k( getNormalisedWaveNumbers() );

#if 0
  // alpha defined as in
  // Bauer, Theory of waves incoherently scattered,
  // Phil. Trans. Roy. Soc. London, 280, 167-191, 1975
  Array1dr alpha( where(k==0.0, 0.0, 1.0/(k*toNormalisedLength(lambdae))) );
  Array1dr alpha2(pow2(alpha));
  Array1dr alpha4(pow2(alpha2));

  // Sanbonmatsu et al., Quantitative comparison of reduced-description
  // particle-in-cell and quasilinear-Zakharov models for parametrically
  // excited Langmuir turbulence, Phys. Plasma, 7, 2824-2841, 2000
  nk2 = where(k==0.0, 1.0/(1.0+Te/Ti),
              (1.0+alpha2)/(1.0+(1.0+Te/Ti)*alpha2));
  // Bauer, Theory of waves incoherently scattered,
  // Phil. Trans. Roy. Soc. London, 280, 167-191, 1975
  nk2 = where(k==0.0, 1.0/(1.0+Te/Ti),
              1.0/(1.0+alpha2)+ alpha4/(1.0+(1.0+Te/Ti)*alpha2)/(1.0+alpha2));
#endif

  // From Guio thesis
  Array1dr a2( pow2(k*toNormalisedLength(lambdae)) );
  Array1dr nk2( 1.0/(1.0+a2)/(1.0+a2+Te/Ti));

  // nk2 = <|\Delta N(k)/N|^2> = 1/N ( correction term)
  nk2 *= pow2(ne)*(1.0/Ne);
  Array1dr nk( sqrt(nk2));
  nk = toNormalisedDensity(nk);

  return nk;
}

Array1dr Mapper::getEkFluctuationsLevel() const
{
  Array1dr k( getNormalisedWaveNumbers() );

#if 0
  // alpha defined as in
  // Bauer, Theory of waves incoherently scattered,
  // Phil. Trans. Roy. Soc. London, 280, 167-191, 1975
  Array1dr alpha( where(k==0.0, 0.0, 1.0/(k*toNormalisedLength(lambdae))) );
  Array1dr alpha2(pow2(alpha));
  Array1dr alpha4(pow2(alpha2));

  // Sanbonmatsu et al., Quantitative comparison of reduced-description
  // particle-in-cell and quasilinear-Zakharov models for parametrically
  // excited Langmuir turbulence, Phys. Plasma, 7, 2824-2841, 2000
  Ek2 = where(k==0.0, 0.0, 1.0/(1+alpha2) +
              alpha2*(1+alpha2)/(1.0+(1.0+Te/Ti)*alpha2)/pow2(1+alpha2));
#endif

  // From Guio thesis
  Array1dr a2( pow2(k*toNormalisedLength(lambdae)) );
#if 0
  Array1dr Ek2( where(k==0.0, 0.0, pow2(lambdae)/2.0));
#elif 0
  Array1dr Ek2( where(k==0.0, 0.0, pow2(lambdae)/(1.0+a2)/2.0));
#else
  Array1dr Ek2( pow2(lambdae)/(1.0+a2)/2.0);
#endif

  Ek2 *= pow2(SI::e/SI::eps0)*pow2(ne)*(1.0/Ne);

  Array1dr Ek( sqrt(Ek2));
  Ek = toNormalisedElectricField(Ek);

#if 0
  dD.resize(N);
  dD = Ek2/(4.0*L_m);
#endif
  return Ek;
}

Array1dr Mapper::getSoundWavesDamping() const
{
  Array1dr k ( getNormalisedWaveNumbers() );

  // Hanssen et al., Numerical test of the weak turbulence approximation to
  // ionospheric Langmuir turbulence, JGR, 97, 12073-12091, 1992, (Eq. 13 p. 12077)
  // Ichimari, Statistical plasma physics, Vol.1, 1992, (Eq. 4.48 p. 118)
  Array1dr wih( sqrt(m_pi/8)*(sqrt(SI::me/mi)+pow(Te/Ti,1.5)*exp(-Te/2/Ti-1.5))*abs(k) );

  // Robinson, Nonlinear wave collapse and strong turbulence,
  // Rev. Mod. Phys., 69, 1997, (Eq. 2.6 p. 511)
  // Be careful typo
  Array1dr wir( sqrt(eta*m_pi/8)*
                (sqrt(SI::me/mi)+eta*pow(.5*Te/Ti,1.5)*exp(-eta*Te/2/Ti))*abs(k) );

  // Numerical method: complex root of dispersion relation
  // (k,w) are in physical units and k's are zero-centered
  Array1dr kp( getPhysicalWaveNumbers() );
  Array1dc w0( zip(kp*cs,
                   -fabs(kp)*sqrt(m_pi/8)*(sqrt(SI::me/mi)+pow(Te/Ti,1.5)*exp(-Te/2/Ti-1.5))*cs,
                   Complex()) );
  Array1dc win(N);
  {
#if 0
    Complex w = getDispersionRelationRoot(kp(0), w0(1));
#else
    Complex w = getDispersionRelationRoot1(kp(0), w0(1));
#endif
    win(0) = Complex( real(w), -imag(w) );
  }
  for (int i=1; i<N/2; ++i) {
    if ( kp(i) ) {
      {
#if 0
        Complex w = getDispersionRelationRoot(kp(i), w0(i));
#else
        Complex w = getDispersionRelationRoot1(kp(i), w0(i));
#endif
        win(i) = Complex( real(w), -imag(w) );
      }
      {
#if 0
        Complex w = getDispersionRelationRoot(kp(N-i), w0(N-i));
#else
        Complex w = getDispersionRelationRoot1(kp(N-i), w0(N-i));
#endif
        win(N-i) = Complex( real(w), -imag(w) );
      }
    } else {
      win(i) = Complex( 0.0, 0.0);
    }
  }
  win = fourier::ifftshift( toNormalisedFreq(win) );
  w0  = fourier::ifftshift( toNormalisedFreq(w0) );

  // sum of collision and Landau damping
#if 1
  Array1dr nui( 0.5*toNormalisedFreq(nuic) + imag(win) );
#else
  Array1dr nui( 0.5*toNormalisedFreq(nuic) + wih );
#endif

  Range isaved(Range(0,adjfn(0,N-1,8),8));
  saveMatlab("nui.m", ios::out|ios::trunc, "k", k(isaved));
  saveMatlab("nui.m", ios::out|ios::app, "wih", wih(isaved));
  saveMatlab("nui.m", ios::out|ios::app, "wir", wir(isaved));
  saveMatlab("nui.m", ios::out|ios::app, "w0", w0(isaved));
  saveMatlab("nui.m", ios::out|ios::app, "win", win(isaved));
  saveMatlab("nui.m", ios::out|ios::app, "nu_i", nui(isaved));

  return nui;
}

Array1dr Mapper::getMaxwellianLangmuirWavesDamping() const
{
  Array1dr k( getNormalisedWaveNumbers() );

  // Hanssen et al., Numerical test of the weak turbulence approximation to
  // ionospheric Langmuir turbulence, JGR, 97, 12073-12091, 1992, (Eq. 14 p. 12077)
  Array1dr wih( where(k==0.0, 0.0, sqrt(m_pi/8)*pow4(1.5)*pow(mi/(eta*SI::me),2.5)/
                      (abs(k)*pow2(k))*exp(-1.125*mi/(eta*SI::me)/pow2(k)-1.5)) );

  // Calculation described in doc/damping.tex
  Array1dr wig( where(k==0.0, 0.0, sqrt(m_pi/8)*pow4(1.5)*pow(mi/(eta*SI::me),2.5)/
                      (abs(k)*pow2(k))*(1.+4/3*pow2(k)*eta*SI::me/mi)*
                      exp(-1.125*mi/(eta*SI::me)/pow2(k)-1.5)) );

  // Numerical method: complex root of dispersion relation
  Array1dr kp( getPhysicalWaveNumbers() );
  Array1dc w0( zip(wpe*(1+1.5*pow2(kp*lambdae)), -fourier::fftshift(toPhysicalFreq(wig)), Complex()) );
  Array1dc win(N);
  {
#if 0
    Complex w = getDispersionRelationRoot(kp(0), w0(0));
#else
    Complex w = getDispersionRelationRoot1(kp(0), w0(0));
#endif
    win(0) = Complex( real(w), -imag(w) );
  }
  for (int i=1; i<N/2; ++i) {
    if ( kp(i) ) {
      {
#if 0
        Complex w = getDispersionRelationRoot(kp(i), w0(i));
#else
        Complex w = getDispersionRelationRoot1(kp(i), w0(i));
#endif
        win(i) = Complex( real(w), -imag(w) );
      }
      {
#if 0
        Complex w = getDispersionRelationRoot(kp(N-i), w0(N-i));
#else
        Complex w = getDispersionRelationRoot1(kp(N-i), w0(N-i));
#endif
        win(N-i) = Complex( real(w), -imag(w) );
      }
    } else {
      win(i) = Complex( real(w0(i)), 0.0);
    }
  }
  win = fourier::ifftshift( toNormalisedFreq(win) );
  w0  = fourier::ifftshift( toNormalisedFreq(w0) );

  // sum of collision and Landau damping
  Array1dr nue( 0.5*toNormalisedFreq(nuec) + wig );

  Range isaved(Range(0,adjfn(0,N-1,8),8));
  saveMatlab("nuet.m", ios::out|ios::trunc, "k", k(isaved));
  saveMatlab("nuet.m", ios::out|ios::app, "wih", wih(isaved));
  saveMatlab("nuet.m", ios::out|ios::app, "wig", wig(isaved));
  saveMatlab("nuet.m", ios::out|ios::app, "w0", w0(isaved));
  saveMatlab("nuet.m", ios::out|ios::app, "win", win(isaved));
  saveMatlab("nuet.m", ios::out|ios::app, "nu_e", nue(isaved));

  return nue;
}

Array1dr Mapper::getLangmuirWavesDamping() const
{
  Array1dr k( getNormalisedWaveNumbers() );

  // Hanssen et al., Numerical test of the weak turbulence approximation to
  // ionospheric Langmuir turbulence, JGR, 97, 12073-12091, 1992
  Array1dr wih( where(k==0.0, 0.0, sqrt(m_pi/8.0)*pow4(3.0/2.0)*pow(mi/(eta*SI::me),2.5)/
                      (abs(k)*pow2(k))*exp(-(9.0/8.0)*mi/(eta*SI::me)/pow2(k)-1.5)) );

  // Calculation described in doc/damping.tex
  Array1dr wig( where(k==0.0, 0.0, sqrt(m_pi/8.0)*pow4(3.0/2.0)*pow(mi/(eta*SI::me),2.5)/
                      (abs(k)*pow2(k))*(1.0+4.0/3.0*pow2(k)*eta*SI::me/mi)*
                      exp(-(9.0/8.0)*mi/(eta*SI::me)/pow2(k)-1.5)) );

  Array1dr alpha(k*2.0/3.0*sqrt((eta*SI::me)/mi));
  Array1dr beta(sqrt(1.0+3.0*pow2(alpha)));
  Array1dr wibf(N), wibg(N);
  wibf = 0.0;
  wibg = 0.0;
  for (unsigned j=0; j<nb_ne.size(); ++j) {
    Real nfrac = nb_ne[j];
    // Forme, Parametric decay of beam-driven Langmuir wave and
    // enhanced ion-acoustic fluctuations, Ann. Geo., 17, 1172-1181, 1999
    Array1dr x( where(k==0.0, 0.0, toNormalisedWaveNumber(kb[j])/k) );
    wibf += where(x==0.0, 0.0, sqrt(m_pi/8.0)*(3.0/2.0)*mi/(eta*SI::me)*abs(x)/x*
                  pow2(x/dub_ub[j])*nfrac*(1.0/dub_ub[j])*(x-1)*
                  exp(-0.5*pow2((1.0-x)/dub_ub[j])));

    // Calculation described in doc/damping.tex
    Real du = dub_ub[j]*fabs(ub[j]);
    Real u = ub[j];
    wibg += where(k==0.0, 0.0, sqrt(m_pi/8.0)*pow3(3.0/2.0)*pow2(mi/(eta*SI::me))/
                  pow2(k)*nfrac*pow2(ve/du)*beta*abs(k)/k*
                  (beta/alpha*(ve/du)-u/du)*
                  exp(-0.5*pow2(beta/alpha*(ve/du)-u/du)));
  }

  // sum of collision and Landau damping
  Array1dr nue( 0.5*toNormalisedFreq(nuec) + wig + wibg );

  Range isaved(Range(0,adjfn(0,N-1,8),8));
  saveMatlab("nue.m", ios::out|ios::trunc, "k", k(isaved));
  saveMatlab("nue.m", ios::out|ios::app, "wih", wih(isaved));
  saveMatlab("nue.m", ios::out|ios::app, "wig", wig(isaved));
  saveMatlab("nue.m", ios::out|ios::app, "wibf", wibf(isaved));
  saveMatlab("nue.m", ios::out|ios::app, "wibg", wibg(isaved));
  saveMatlab("nue.m", ios::out|ios::app, "nu_e", nue(isaved));

  return nue;
}

Complex Mapper::getPump(Real t) const
{
  Real phase = dOmega_n*t;
#if 0
  cout << "t=" << t << " phase=" << phase
       << " cos(phase)=" << cos(phase) << " sin(phase)=" << sin(phase) << endl;
#endif
  return Complex(Epump_n*cos(phase), -Epump_n*sin(phase));
}

Real Mapper::getN2Wconstant() const
{
  if (physical_units) {
    // Parseval'theorem  \sum |E_j|^2 = 1/nfft \sum |E_k|^2
    // N is defined as Eq. 34 in Payne et al, 1983
    // N = \sum |E_j/\epsilon|^2 (L_m/\chi)/nfft
    //   = \sum |E_k/\epsilon|^2 (L_m/\chi)/nfft^2
    // E_rms^2 = 1/nfft \sum |E_j|^2 = N\epsilon^2\chi/L_m
    // Eq. 2.21 from Robinson 1997
    // W = \epsilon_0 |E_rms|^2/(2n_e K_b T_e)
    //   = \epsilon_0 N\epsilon^2\chi/L_m/(2n_e K_b T_e)
    // K = \epsilon_0 /(2n_e K_b T_e) \epsilon^2\chi / L_m
    // K = \epsilon_0 /(2n_e K_b T_e) \epsilon^2 / L
    // K = 8/3 \eta^2 (me/mi) / L
    // K = 8/3 (ce+ci*Ti/Te)^2 (me/mi) / L
#if 0
    std::cout << "N2W=" <<  8.0/3.0*pow2(eta)*(SI::me/mi)/L << " , "
              << SI::eps0/(2.0*ne*SI::Kb*Te)*pow2(toPhysicalElectricField(1.0))/L
              << std::endl;
#endif
    return SI::eps0/(2.0*ne*SI::Kb*Te)*pow2(toPhysicalElectricField(1.0))/L;
  } else {
    // Assuming Te>>Ti and electron-to-proton mass ratio (\sim 1/1836)
    return 8.0/3.0*(SI::me/SI::mp)/L;
  }
}

Complex Mapper::getDispersionRelation(Real k, Complex w) const
{
  Real li2 = SI::Kb*SI::eps0/pow2(SI::e)*Ti/ne;
  Complex ei = W( w / fabs(k*vi) ) / (pow2(k)*li2);

  Real le2 = SI::Kb*SI::eps0/pow2(SI::e)*Te/ne;
  Complex ee = W( w / fabs(k*ve) ) / (pow2(k)*le2);

  return ei + ee + 1.0;
}


Complex Mapper::getDispersionRelation1stDerivative(Real k, Complex w) const
{
#if 0
  const Real dz = 1e-6;
  return (getDispersionRelation(k, w*(1.+dz))-getDispersionRelation(k, w*(1.-dz)))/(2.*dz*w);
#else

  Real li2 = SI::Kb*SI::eps0/pow2(SI::e)*Ti/ne;
  Complex deidw = dWdz( w / fabs(k*vi) ) / fabs(k*vi) / (pow2(k)*li2);

  Real le2 = SI::Kb*SI::eps0/pow2(SI::e)*Te/ne;
  Complex deedw = dWdz( w / fabs(k*ve) ) / fabs(k*ve) / (pow2(k)*le2);

  return deidw + deedw;
#endif
}

Complex Mapper::getDispersionRelation2ndDerivative(Real k, Complex w) const
{
  Real li2 = SI::Kb*SI::eps0/pow2(SI::e)*Ti/ne;
  Complex d2eidw2 = d2Wdz2( w / fabs(k*vi) ) / pow2(k*vi) / (pow2(k)*li2);

  Real le2 = SI::Kb*SI::eps0/pow2(SI::e)*Te/ne;
  Complex d2eedw2 = d2Wdz2( w / fabs(k*ve) ) / pow2(k*ve) / (pow2(k)*le2);

  return d2eidw2 + d2eedw2;
}



Complex Mapper::getDispersionRelationRoot(Real k, Complex wi) const
{
  const Real epsilon = 1e-6;
  const int maxIter = 2000;
  Complex w(wi);
  Complex F( getDispersionRelation(k, w) );
  Complex dFdw( getDispersionRelation1stDerivative(k, w) );
  Complex d2Fdw2( getDispersionRelation2ndDerivative(k, w) );
  Complex dw;

  int i = 0;
#if 0
  cout << "k=" << k << " w=" << w << " W(w)=" << F << endl;
#endif
  do {
#if 0
    // Newton-Raphson to get the root of the dispersion relation
    dw = F / dFdw;
#else
    // Halley's iteration formula for f(z)=0
    // Description at http://mathworld.wolfram.com/HalleysMethod.html
    // x_{n+1} = x_n - 2 f(x_n) f'(x_n) / (2 f'(x_n)^2-f(x_n)f''(x_n) )
    //dw =  2.0 * F * dFdw / (2.0 * pow2(dFdw) - F * d2Fdw2);
    dw =  F / (dFdw - 0.5 * F * d2Fdw2 / dFdw);
#endif
    w -= dw;
    if (imag(w) > 0.0) {
      cout << "WARNING Imag(w)>0: (k,w) = "
           << std::showpos <<  k << ", " << w << std::noshowpos << endl;
      w = conj(w);
    }
    F = getDispersionRelation(k, w);
    dFdw = getDispersionRelation1stDerivative(k, w);
    d2Fdw2 = getDispersionRelation2ndDerivative(k, w);
    ++i;
  } while (abs(dw) > epsilon * abs(w) && abs(F) > epsilon && i < maxIter);

  if (i == maxIter) {
    ostringstream os;
    os << "No convergence achieved after " << i << " iterations\n"
       << "Initial values k, w = " << k << ", " << wi << ", " << getDispersionRelation(k, wi) << "\n"
       << "Last value w, dw, disp. = " << w << ", " << dw << ", " << F;
    throw ClassException("Mapper", os.str());
  }

#if 0
  cout << "0 iter=" << i << " k, w=" << k << "," << w <<
       " W(w)=" << F << endl;
#endif
  return w;
}

Complex Mapper::getDispersionRelationRoot1(Real k, Complex wi) const
{
  const Real epsilon = 1e-6;
  const Real h = 1e-5;
  const int maxIter = 2000;
  Complex dw;
  Complex w[3] = { wi, wi*(1.0+h), wi*(1.0+h)*(1+h) };
  Complex F[3] = { getDispersionRelation(k, w[0]),
                   getDispersionRelation(k, w[1]),
                   getDispersionRelation(k, w[2])
                 };
  int i = 0;
#if 0
  cout << "k=" << k << " w=" << w[0] << " W(w)=" << F[0] << endl;
#endif
  do {
#if 0
    // Newton method without derivative
    Complex D ( (F[0]-F[1])/(w[0]-w[1]) );
    dw = F[0] / D;
    w[1] = w[0];
    w[0] -= dw;
    F[1] = F[0];
    F[0] = getDispersionRelation(k, w[0]);
#else
    // Halley's iteration method without derivative
    // Eq. 23-27 in Gillan et al. Comp. Phys. Comm. 175 (2006) 304-313
    Complex dF[2] = { (F[0]-F[1])/(w[0]-w[1]), (F[1]-F[2])/(w[1]-w[2]) };
    Complex D ( dF[0] - F[1] * (dF[0] - dF[1]) / (w[0] - w[2]) / dF[0] );
    dw = F[0] / D;
    w[2] = w[1];
    w[1] = w[0];
    w[0] -= dw;
    F[2] = F[1];
    F[1] = F[0];
    F[0] = getDispersionRelation(k, w[0]);
#endif
    if (imag(w[0]) > 0.0) {
      cout << "WARNING Imag(w)>0: (k,w) = "
           << std::showpos << k << ", " << w[0] << std::noshowpos << endl;
      w[0] = conj(w[0]);
    }
    ++i;
  } while (abs(dw) > epsilon * abs(w[0]) && abs(F[0]) > epsilon && i<maxIter);

  if (i == maxIter) {
    ostringstream os;
    os << "No convergence achieved after " << i << " iterations\n"
       << "Initial values k, w, disp. = " << k << ", " << w << ", " << getDispersionRelation(k, w[0]) << "\n"
       << "Last value w, dw, disp. = " << w[0] << ", " << dw << ", " << F[0];
    throw ClassException("Mapper", os.str());
  }

#if 0
  cout << "1 iter=" << i << " k, w=" << k << "," << w[0]
       << " W(w)=" << F[0] << endl;
#endif
  return w[0];
}

Array1dr Mapper::initDiffusion()
{
  Array1dr k( getPhysicalWaveNumbers() );
  vphi.resize(N);
  // grid of phase velocities
  Array1dr alpha(k*lambdae);
  Array1dr wr(wpe*sqrt(1.0+3.0*pow2(alpha)));
  vphi = where(k==0.0, 0.0, wr/k);

  Fe.resize(N);
  // Maxwellian background
  Fe = where(k==0.0, 0.0, 1.0/sqrt(m_2pi)/ve*exp(-0.5/pow2(alpha)-1.5));

#if SEP // Fem separated
  Fem.resize(N);
  Fem = Fe;
  Fe  = 0.0;
#endif

  // Beams
  Array1dr beta(sqrt(1.0+3.0*pow2(alpha)));
  for (unsigned j=0; j<nb_ne.size(); ++j) {
    Real nfrac = nb_ne[j];
    Real du = dub_ub[j]*fabs(ub[j]);
    Real u = ub[j];
    Fe += where(k==0.0, 0.0, nfrac/sqrt(m_2pi)/du*
                exp(-0.5*pow2(beta/alpha*(ve/du)-u/du)));
  }

  Fet0.resize(N);
  Fet0 = Fe;


  dkdvphi.resize(N);
  dkdvphi = where(k==0.0, 0.0, -k*abs(k)/pow2(wpe)*wr);

  dk = (k(N-1)-k(0))/(N-1);

  dFe.resize(N);
  gamma.resize(N);

  dampingCoeff.resize(N);
  // Initialise damping coefficient
  dampingCoeff = where(k==0.0, 0.0, -m_pi/2.0*pow2(wpe/k)*wr);
  dampingCoeff = toNormalisedFreq(dampingCoeff);

#if SEP // Fem separated
  Fe = Fet0+Fem;
#endif

  Range rp(0,N/2);
  deriveFeWithVphi(rp);
  Range rm(N/2+1,N-1);
  deriveFeWithVphi(rm);

#if SEP // Fem separated
  Fe = Fet0;
#endif

  gamma = dampingCoeff*dFe;
  Array1dr nue( fourier::ifftshift( Array1dr(0.5*toNormalisedFreq(nuec) + gamma) ));

  // Initialise diffusion constant
  diffusionCoeff.resize(N);
#if 0
  diffusionCoeff = pow2(m_2pi)*pow2(SI::e/SI::me)*abs(k/wpe)/L_m;
#endif

  // Eq. 10 from Sanbonmatsu et al., Quantitative comparison of reduced-
  // description particle-in-cell and quasilinear-Zakharov models for
  // parametrically excited Langmuir turbulence,
  // Phys. Plasma, 7, 2824-2841, 2000
  diffusionCoeff = where(k==0.0, 0.0, pow2(m_2pi*SI::e/SI::me)/abs(vphi)/L_m);
#if 0
  diffusionCoeff = where(k==0.0, 0.0, pow2(m_2pi*SI::e/SI::me)/abs(vphi));
  diffusionCoeff = where(k==0.0, 0.0, pow2(m_2pi*SI::e/SI::me)/abs(vphi)/L_m);
#endif

  D.resize(N);

#if 0 // buffer to estimate time mean value of <E(k)^2>
  Ekt2.resize(N,nt);
  Ekt2 = 0.0;
#endif

#if 0
  dD = dD*diffusionCoeff;
#endif

  // Initialise ranges to integrate
  setIntegrationBoundaries();

#if 0 // alternative diffusion term
  wk.resize(N);
  wk = wpe*(1+1.5*abs(k)*lambdae);
  wk = where(k < 0.0, -wk, wk);
  gk.resize(N);
  gk = toPhysicalFreq(gamma);
  gk = where(k < 0.0, -gk, gk);
  gk2.resize(N);
  gk2 = gk*gk;
  for (int i=ipRange.first(); i<=ipRange.last(); ++i) {
    Complex w = Complex(wk(i), gk(i));
    cout << "k, w = " << k(i) << ", (" << wk(i) << "," << gk(i) << ") ";
    w = getDispersionRelationRoot(k(i), w);
    wk(i) = real(w);
    gk(i) = imag(w);
    cout << "w = (" << wk(i) << "," << gk(i) << ")" << endl;
  }

  for (int i=inRange.first(); i<=inRange.last(); ++i) {
    Complex w = Complex(wk(i), gk(i));
    cout << "w = " << wk(i) << "," << gk(i) << " ";
    w = getDispersionRelationRoot(k(i), w);
    wk(i) = real(w);
    gk(i) = imag(w);
    cout << "w = " << wk(i) << "," << gk(i) << endl;
  }
#endif

  Range isaved(Range(0,adjfn(0,N-1,8),8));
  saveMatlab("diffusion.m", ios::out|ios::trunc, "kd", k(isaved));
  saveMatlab("diffusion.m", ios::out|ios::app, "Fe", Fe(isaved));
  saveMatlab("diffusion.m", ios::out|ios::app, "dFe", dFe(isaved));
  saveMatlab("diffusion.m", ios::out|ios::app, "nue", nue(isaved));

  return nue;
}

void Mapper::updateDiffusion(const Array1dc &Ek, Array1dr &nue)
{
  // Calculate D(v,t) from Eq. 10 of
  // Sanbonmatsu et al., Quantitative comparison of reduced-description
  // particle-in-cell and quasilinear-Zakharov models for parametrically
  // excited Langmuir turbulence, Phys. Plasma, 7, 2824-2841, 2000

  Array1dr Ek2(pow2(real(Ek))+pow2(imag(Ek)));
  Ek2 = Ek2*1.0/pow2(N);
  Ek2 = fourier::fftshift( toPhysicalElectricField2(Ek2) );

#if 0 // buffer to estimate time mean value of <E(k)^2>
  Ekt2(Range::all(), it) = Ek2;
  ++it;
  if (!bufferComplete && it == nt)
    bufferComplete = true;
  it %= nt;
  Ek2 = blitz::sum(Ekt2, blitz::secondIndex);
  double s = (bufferComplete ? 1.0/nt : 1.0/it);
#if 0
  cout << "it=" << it << "s=" << 1.0/s << endl;
#endif
  Ek2 *= s;
#endif

#if 0
  Ek2 = fourier::binom8filter(Ek2, Real(0));
#endif

  D = diffusionCoeff*Ek2;

  Array1dr k ( getPhysicalWaveNumbers() );
  Range isaved(Range(0,adjfn(0,N-1,8),8));
  saveMatlab("D.m", ios::out|ios::trunc, "k", k(isaved));
  saveMatlab("D.m", ios::out|ios::app, "Di", D(isaved));

#if 0 // alternative diffusion term
  Array1dr D1;
  D1.resize(N);
  D1 = 0.0;
  for (int i=ipRange.first(); i<=ipRange.last(); ++i) {
    D1(i) = 2*pow2(SI::e/SI::me)*sum(pow2(k(ipRange))*Ek2(ipRange)*abs(gk(ipRange))/
                                     (pow2(k(ipRange)*vphi(i)-wk(ipRange))+gk2(ipRange)));
  }
  for (int i=inRange.first(); i<=inRange.last(); ++i) {
    D1(i) = 2*pow2(SI::e/SI::me)*sum(pow2(k)*Ek2*abs(gk)/(pow2(k*vphi(i)-wk)+gk2));
  }
#endif

#if SEP // Fem separated
  Fe = Fet0;
  Fet0 = Fem+Fet0;
#endif
  // integrate diffusion equation
  crankNicholson(inRange);
  crankNicholson(ipRange);
#if SEP // Fem separated
  Fet0 = Fe;
  Fe = Fet0+Fem;
#endif

#if 0
  ftcsLax(inRange);
  ftcsLax(ipRange);
#endif

  deriveFeWithVphi(in1Range);
  deriveFeWithVphi(ip1Range);

  gamma = dampingCoeff*dFe;
  nue = fourier::ifftshift( Array1dr(0.5*toNormalisedFreq(nuec) + gamma) );

#if defined(HAVE_GNUPLOT)
  gnuplot();
#endif
}

void Mapper::setIntegrationBoundaries()
{
  using blitz::tensor::i;

  Array1di l(vphi.size());
  Array1di lmax(vphi.size());
  Array1di lmin(vphi.size());

  Real vmin = (vphi(N/2+M) > 1.3*ve ? vphi(N/2+M) : 1.3*ve);
  Real vmax = 333.0*ve;

  l = where(vphi > vmin && vphi < vmax, 1, 0);
  lmax = l*i;
  lmin = l*(N-1-i);
  ipRange.setRange(N-1-max(lmin), max(lmax));
  ip1Range.setRange(N-1-max(lmin)+1, max(lmax)-1);

  l = where(vphi < -vmin && vphi > -vmax, 1, 0);
  lmax = l*i;
  lmin = l*(N-1-i);
  inRange.setRange(N-1-max(lmin), max(lmax));
  in1Range.setRange(N-1-max(lmin)+1, max(lmax)-1);

  cout << "ipRange = " << ipRange << endl;
  cout << "inRange = " << inRange << endl;
}

void Mapper::deriveWithK(const Array1dr &y, Array1dr &dy)
{
  dy(0) = (y(1)-y(0))/dk;

  int n = y.rows();
  Range i(1,n-2);
  dy(i) = 0.5*(y(i+1)-y(i-1))/dk;

  dy(n-1) = (y(n-1)-y(n-2))/dk;
}

void Mapper::deriveFeWithVphi(Range &range)
{
  Array1dr yr(Fe(range));
  Array1dr dyr(dFe(range));
  deriveWithK(yr, dyr);
  dyr = dyr*dkdvphi(range);
}


void Mapper::ftcsLax(Array1dr &u, Array1dr &d, Array1dr &dkdv)
{
  // FTCS (Forward-Time Central-Space) scheme corrected by Lax method
  // for advection equation

  // Courant condition for stability
  Array1dr c(abs(dkdv)*d*h/pow2(dk));
  if (max(c)>1)
    throw ClassException("Mapper",
                         "Courant condition not fullfilled in ftcsLax");

  int n = u.rows();
  Range i(1,n-2);

  Array1dr u1(n);
  u1 = u;
  u(i) = 0.5*(u1(i-1)+u1(i+1))+
         0.5*dkdv(i)*h/dk*(d(i+1)*u1(i+1)-d(i-1)*u1(i-1));
}

void Mapper::ftcsLax(Range &range)
{
  Array1dr fe(Fe(range));
  Array1dr d(vphi(range)*dD(range)/ve2);
  Array1dr dkdv(dkdvphi(range));
  ftcsLax(fe, d, dkdv);
}

void Mapper::crankNicholson(Array1dr &u, Array1dr &d, Array1dr &dkdv)
{
  // Crank-Nicholson scheme (semi-implicit method)
  // with Dirichlet boundary condition
  int n = u.rows();
  Array1dr di(n-1);
  Array1dr dkdvi(n-1);
  Array1dr a(n), b(n), c(n);
  Array1dr u1(n);
  Real alpha = h/(2.0*dk*dk);

  Range i(0,n-2);
  di = 0.5*(d(i)+d(i+1)); // D_{j-1/2}^n
  dkdvi = 0.5*(dkdv(i)+dkdv(i+1)); // d k_{j+1/2) / dv
  di = di*dkdvi;

  i.setRange(0,n-3);

  a(0) = a(n-1) = 0.0;
  a(i+1) = -dkdv(i+1)*alpha*di(i);

  b(0) = b(n-1) = 1.0;
  b(i+1) = 1.0+dkdv(i+1)*alpha*(di(i)+di(i+1));

  c(0) = c(n-1) = 0.0;
  c(i+1) = -dkdv(i+1)*alpha*di(i+1);

  u1(0) = u(0);
  u1(n-1) = u(n-1);
  u1(i+1) = (dkdv(i+1)*alpha*di(i)*u(i)+
             (1.0-dkdv(i+1)*alpha*(di(i)+di(i+1)))*u(i+1)+
             dkdv(i+1)*alpha*di(i+1)*u(i+2));

  tridiag(a, b, c, u1, u);
}

void Mapper::crankNicholson(Array1dr &u, Array1dr &d, Array1dr &dkdv,
                            Array1dr &v, Array1dr &fet0)
{
  // Crank-Nicholson scheme (semi-implicit method)
  // with Dirichlet boundary condition
  // and advective term
  int n = u.rows();
  Array1dr di(n-1);
  Array1dr dkdvi(n-1);
  Array1dr a(n), b(n), c(n);
  Array1dr u1(n);
  Real alpha = h/(2.0*dk*dk);
  Real advectConst = 2.0;

  Range i(0,n-2);
  di = 0.5*(d(i)+d(i+1));
  dkdvi = 0.5*(dkdv(i)+dkdv(i+1));
  di = di*dkdvi;

  i.setRange(0,n-3);

  a(0) = a(n-1) = 0.0;
  a(i+1) = -dkdv(i+1)*alpha*di(i);

  b(0) = b(n-1) = 1.0;
#if 1
  b(i+1) = 1.0+dkdv(i+1)*alpha*(di(i)+di(i+1))+abs(v(i+1))*h/(2.0*advectConst*L_m);
#else
  b(i+1) = 1.0+dkdv(i+1)*alpha*(di(i)+di(i+1));
#endif

  c(0) = c(n-1) = 0.0;
  c(i+1) = -dkdv(i+1)*alpha*di(i+1);

  u1(0) = u(0);
  u1(n-1) = u(n-1);
#if 1
  u1(i+1) = (dkdv(i+1)*alpha*di(i)*u(i)+
             (1.0-dkdv(i+1)*alpha*(di(i)+di(i+1))-
              abs(v(i+1))*h/(2.0*advectConst*L_m))*u(i+1)+
             dkdv(i+1)*alpha*di(i+1)*u(i+2))+abs(v(i+1))*h/(advectConst*L_m)*fet0(i+1);
#else
  u1(i+1) = (dkdv(i+1)*alpha*di(i)*u(i)+
             (1.0-dkdv(i+1)*alpha*(di(i)+di(i+1)))*u(i+1)+
             dkdv(i+1)*alpha*di(i+1)*u(i+2))-(v(i+1))*h/(advectConst*L_m)*(u(i+1)-fet0(i+1));
#endif

  tridiag(a, b, c, u1, u);
}

void Mapper::crankNicholson(Range &range)
{
  Array1dr fe(Fe(range));
  Array1dr d(D(range));
  Array1dr dkdv(dkdvphi(range));
#if 1 // without advection term
  crankNicholson(fe, d, dkdv);
#else // with advection term
  Array1dr v(vphi(range));
  Array1dr fet0(Fet0(range));
  crankNicholson(fe, d, dkdv, v,fet0);
#endif
}


void Mapper::tridiag(Array1dr &a, Array1dr &b, Array1dr &c,
                     Array1dr &r, Array1dr &u)
{
  // Solver for a vector u, the tridiagonal linear set given by
  // A u = r where A is tridiagonal a2...an, b1...bn, c1...cn-1
  Real bet;
  int n = r.rows();
  Array1dr gam(n); // Workspace vector

  if (b(0) == 0.0)
    throw ClassException("Mapper", "Error in tridag");

  u(0) = r(0)/(bet=b(0));
  for (int j=1; j<=n-1; ++j) {
    // Decomposition and forward substitution
    gam(j) = c(j-1)/bet;
    bet = b(j)-a(j)*gam(j);
    if (bet == 0.0)
      throw ClassException("Mapper", "Error in tridag");
    u(j) = (r(j)-a(j)*u(j-1))/bet;
  }
  for (int j=n-2; j>=0; j--) {
    // Backsubstitution
    u(j) -= gam(j+1)*u(j+1);
  }
}

Complex Mapper::wofz(const Complex z) const
{
  // complex scaled complementary error function (FADDEEVA function)

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
  if (xabs>rmaxreal || yabs>rmaxreal)
    throw ClassException("Mapper",
                         "Error in wofz\n\txabs>rmaxreal || yabs>rmaxreal");

  Real qrho = x*x + y*y;

  Real xabsq = xabs*xabs;
  Real xquad = xabsq - yabs*yabs;
  Real yquad = 2.0*xabs*yabs;

  Real u2 = 0;
  Real v2 = 0;

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
    Real h2 = 0;
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

    Real qlambda = 0;
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
      if (yquad>rmaxgoni || xquad>rmaxexp)
        throw ClassException("Mapper",
                             "Error in wofz\n\tyquad>rmaxgoni || xquad>rmaxexp");

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

  return Complex(u,v);
}

Complex Mapper::W(const Complex z) const
{
  // Plasma Dispersion Function W defined by
  //
  //                     / inf
  //            1        |     x exp(-x^2/2)     \__/
  // W(y)= -----------   |    -------------- dx,  \/ y / Im(y) > 0
  //       (2 pi)^1/2    |        x - y
  //                    / -inf
  //
  // related to the Z (Fried-Conte) function
  //              y       y          1       y
  // W(y) = 1 + ----- Z(------) = - --- Z'(-----)
  //            2^1/2   2^1/2        2     2^1/2
  // or to the wofz (Faddeeva) function
  //                           y          y
  // W(y) = 1 + i \sqrt{\pi} ----- wofz(------)
  //                         2^1/2      2^1/2
  //
  // Z(z) = i \sqrt{\pi} wofz(z)
  Complex z2 = m_sqrt1_2*z;
  return Complex(1.0, 0.0) + Complex(0.0, m_sqrt_pi)*z2*wofz(z2);
}

Complex Mapper::dWdz(const Complex z) const
{
  // dW/dz = (1/z - z) W(z) - 1/z
  // dW/dz = 1/sqrt{2} Z(z/sqrt{2})-z -z^2 Z(z/sqrt{2})
  // dW/dz = i\sqrt{\pi/2} wofz(z/sqrt{2})-z -i z^2 \sqrt{\pi} wofz(z/sqrt{2})
#if 0
  {
    const Real h=1e-5;
    Complex Z (Complex(0.0, m_sqrt_pi) * wofz(m_sqrt1_2*z));
    cout << "z=" << z << ", dWdz=" << (W(z+h)-W(z))/h << ", " <<
         ( 1.0 - pow2(z) ) * m_sqrt1_2 * Z - z << endl;
  }
#endif
#if 1
  Complex Z (Complex(0.0, m_sqrt_pi) * wofz(m_sqrt1_2*z));
  return  ( 1.0 - pow2(z) ) * m_sqrt1_2 * Z - z;
#else
  Complex zinv = 1.0 / z;
  return (zinv -z)*W(z) - zinv;
#endif
}

Complex Mapper::d2Wdz2(const Complex z) const
{
  // d^2W/dz^2 = -2 +z^2 +(z^3 -3z)/\sqrt{2} Z(z/\sqrt{2})
  // Z(z) = i \sqrt{\pi} wofz(z)
#if 0
  {
    const Real h=1e-5;
    Complex Z (Complex(0.0, m_sqrt_pi) * wofz(m_sqrt1_2*z));
    cout << "z=" << z << ", d2Wdz2=" << (dWdz(z+h)-dWdz(z))/h << ", " <<
         -2.0 + pow2(z) + (pow2(z) - 3.0) * z * m_sqrt1_2 * Z << endl;
  }
#endif
  Complex Z (Complex(0.0, m_sqrt_pi)*wofz(m_sqrt1_2*z));
  return -2.0 + pow2(z) + (pow2(z) - 3.0) * z * m_sqrt1_2 * Z;
}

#if defined(HAVE_GNUPLOT)
void Mapper::gnuplot()
{
  static gnuplot::Gnuplot g1("lines");
  //	static gnuplot::Gnuplot g2("lines");
  //	static gnuplot::Gnuplot g3("lines");
  static int i = 0, istart=0, imodulo=1000;
  if (i >= istart && i % imodulo == 0) {
    //	g1.reset();
    {
      Array1dr x(vphi(ipRange)/ve);
      Array1dr y(Fe(ipRange));

      y = Fe(ipRange);
      if (i==istart) {
        g1.sendCmd("set xrange [4.5:30]");
        g1.sendCmd("set grid");
        g1.setXLabel("v/v_e");
        g1.setTitle("F_e");
      }
      g1.plot(x, y);
#if 0

      y = abs(dFe(ipRange));
      if (i==istart) {
        g2.sendCmd("set xrange [4.5:30]");
        g2.sendCmd("set grid");
        g2.setXLabel("v/v_e");
        g2.setTitle("dF_e");
      }
      g2.plot(x, y);
#endif
#if 0

      y = D(ipRange);
      if (i==istart) {
        g3.sendCmd("set xrange [4.5:30]");
        g3.sendCmd("set grid");
        g3.setXLabel("v/v_e");
        g3.setTitle("D");
      }
      g3.plot(x, y);
#endif

    }

#if 0
    {
      Array1dr x(abs(vphi(inRange))/ve);
      Array1dr y(Fe(inRange));
      g1.sendCmd("set xrange [4.5:30]");
      g1.plot(x,y,"Fe-");
    }
#endif

  }
  ++i;
}
#endif
