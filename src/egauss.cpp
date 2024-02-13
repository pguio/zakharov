/**************************************************************************
 *
 * $Id: egauss.cpp,v 1.65 2011/03/26 07:47:28 patrick Exp $
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

#include <modes-factory.h>
#include <egauss.h>

using blitz::exp;
using blitz::imag;
using blitz::min_exponent10;
using blitz::pow2;
using blitz::real;
using blitz::where;

using std::cout;
using std::endl;
using std::ostream;

namespace {
  RegisterInFactory<InitMode, EGauss> registerMe("egauss");
}

namespace factory {
  void dummyEgauss()
  {}
}



ostream& operator<<(ostream& os, const EGauss &e)
{
  e.printOn(os);
  return os;
}


EGauss::EGauss(int nargs, char *args[]) : InitMode(nargs, args, "EGauss"),
  E0(0.0), Erms(0.0), nrms(0.0),
  k0(8), W(0.04), Delta(1.0)
{
  initParsing(nargs, args);
  parseParam();
}

EGauss::~EGauss()
{}


void EGauss::initialise(Mapper &mapper)
{
  InitMode::initialise(mapper);

  source_nk = nrms*N;
  source_nk(0) = 0.0;

  source_Ek = Erms*N;
  source_Ek(0) = E0*N;
}

void EGauss::initEndn(Mapper &mapper, FT1d &fft,
                      Array1dc &Ej, Array1dc &nj, Array1dc &dnj,
                      Array1dc &Ek, Array1dc &nk, Array1dc &dnk)
{
  InitMode::initEndn(mapper, fft, Ej, nj, dnj, Ek, nk, dnk);

  Array1dc z1(N);
  Array1dc gauss(N);
  Array1dr k( mapper.getNormalisedWaveNumbers() );

  // Gaussian electric field
  Real E_0 = sqrt(2.0*W/Delta/sqrt(m_2pi)/m_2pi);
#if 0
  cout << "E_0=" << E_0 << endl;
#endif

  for (int i=1; i<=M; ++i) {
    z1(i)   = Complex(0.0, m_2pi*rngE.random());
    z1(N-i) = Complex(0.0, m_2pi*rngE.random());
  }
  gauss = E_0*exp(-pow2((k-k0)/Delta))*complexExp(z1);

  Real z = 1.0;
  Real dmin = pow(10.0, (Real)min_exponent10(z));
#if 0
  cout << "dmin = " << dmin << endl;
#endif
  gauss = where(pow2(real(gauss))+pow2(imag(gauss)) <= dmin, 0.0, gauss);
  gauss = gauss*(Real)N;

  Ek += gauss;
  fft.inverse(Ek, Ej);
}

void EGauss::printOn(ostream& os) const
{
  InitMode::printOn(os);
  os
      << '\n'
      << "E0    = " << E0 << '\n'
      << "Erms  = " << Erms << '\n'
      << "nrms  = " << nrms << '\n'
      << "k0    = " << k0 << '\n'
      << "W     = " << W << '\n'
      << "Delta = " << Delta;
}


void EGauss::initParsing(int nargs, char *args[])
{
  const string id(PACKAGE);
  const string version(
    VERSION "\n"
    "$Id: egauss.cpp,v 1.65 2011/03/26 07:47:28 patrick Exp $"
    "\n");
  const string copyright(ZAK_COPYRIGHT);

  registerClass("EGauss");
  registerPackage(id, version, copyright);

  using parser::types::real;

  insertOption(_E0   , "E0"   , real, "Parameter E0 "          , Any(E0));
  insertOption(_Erms , "Erms" , real, "Parameter Erms"         , Any(Erms));
  insertOption(_nrms , "nrms" , real, "Parameter nrms"         , Any(nrms));

  insertOption(_k0   , "k0"   , real, "Central wave vector"    , Any(k0));
  insertOption(_W    , "W"    , real, "Langmuir wave intensity", Any(W));
  insertOption(_Delta, "Delta", real, "Wave vector width"      , Any(Delta));
}

void EGauss::parseParam()
{
  parseOption(_E0, E0);
  parseOption(_Erms, Erms);
  parseOption(_nrms, nrms);

  parseOption(_k0, k0);
  parseOption(_W, W);
  parseOption(_Delta, Delta);
}
