/**************************************************************************
 *
 * $Id: default-mode.cpp,v 1.32 2020/04/22 13:34:53 patrick Exp $
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

#include <default-mode.h>

using std::ostream;

namespace factory {
  void dummyDefault()
  {}
}


ostream& operator<<(ostream& os, const InitMode &mode)
{
  mode.printOn(os);
  return os;
}


InitMode::InitMode(int nargs, char *args[], const string modename) :
  Parser(nargs, args), ModeName(modename), mseed(1)
{
  initParsing(nargs, args);
  parseParam();

  rngE.seed(mseed(0));
  rngn.seed(mseed(1));
}


InitMode::~InitMode()
{}


void InitMode::initialise(Mapper &mapper)
{
  N = mapper.getGridSize();
  M = mapper.getZeroPaddingSize();

  source_nk.resize(N);
  source_nk = Complex(0.0,0.0);

  source_Ek.resize(N);
  source_Ek = Complex(0.0,0.0);
}

void InitMode::initEndn(Mapper &mapper, FT1d &fft,
                        Array1dc &Ej, Array1dc &nj, Array1dc &dnj,
                        Array1dc &Ek, Array1dc &nk, Array1dc &dnk)
{
  // source_Ek and source_nk should have been initialised

  Array1dc z1(N);

  // Electric field thermal fluctuation
  Ek.resize(N);
  z1 = Complex(0.0, 0.0);
  for (int i=1; i<=M; ++i) {
    z1(i)   = Complex(0.0, m_2pi*rngE.random());
    z1(N-i) = Complex(0.0, m_2pi*rngE.random());
  }
  Ek = complexExp(z1)*source_Ek;

  Ej.resize(N);
  fft.inverse(Ek, Ej);

  // Density thermal fluctuation
  nk.resize(N);
  z1 = Complex(0.0, 0.0);
  // nk( k=0 ) = 0
  for (int i=1; i<=M; ++i) {
    // Make the density spectral fluctuation such that
    // nk(-k) = conj(nk(k))
    z1(i)   = Complex(0.0, m_2pi*rngn.random());
    z1(N-i) = -z1(i);
  }
  nk = complexExp(z1)*source_nk;

  nj.resize(N);
  fft.inverse(nk, nj);

  // Density derivative at t=0
  dnk.resize(N);
  dnk = Complex(0.0,0.0);
  dnj.resize(N);
  fft.inverse(dnk, dnj);
}

void InitMode::getSources(Array1dc &delta_nk, Array1dc &delta_Ek)
{
  delta_nk.resize(N);
  delta_nk = source_nk;
  delta_Ek.resize(N);
  delta_Ek = source_Ek;
}
void InitMode::printOn(ostream& os) const
{
  os << parser::header("Mode setup")
     << "Mode  = " << ModeName << '\n'
     << "mseed = " << mseed;
}

void InitMode::initParsing(int nargs, char *args[])
{
  const string id(PACKAGE);
  const string version(
    VERSION "\n"
    "$Id: default-mode.cpp,v 1.32 2020/04/22 13:34:53 patrick Exp $"
    "\n");
  const string copyright(ZAK_COPYRIGHT);

  registerClass("DefaultMode");
  registerPackage(id, version, copyright);

  using parser::types::intVect;

  insertOption(_mseed, "mseed", intVect,
               "Seeds for the initialisation mode RNGs (E,n)", Any(mseed));
}

void InitMode::parseParam()
{
  parseOption(_mseed, mseed);
}

