/**************************************************************************
 *
 * $Id: solitons.cpp,v 1.48 2011/03/26 07:47:28 patrick Exp $
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
#include <solitons.h>

using blitz::cosh;
using blitz::floor;
using blitz::pow2;
using blitz::sinh;
using blitz::zip;

using std::ios;
using std::ostream;


namespace {
  RegisterInFactory<InitMode, Solitons> registerMe("solitons");
}

namespace factory {
  void dummySolitons()
  {}
}

ostream& operator<<(ostream& os, const Solitons &s)
{
  s.printOn(os);
  return os;
}

Solitons::Solitons(int nargs, char *args[]) : InitMode(nargs, args, "Solitons")
{
  initParsing(nargs, args);
  parseParam();
  checkParam();
}

Solitons::~Solitons()
{}


void Solitons::initialise(Mapper &mapper)
{
  InitMode::initialise(mapper);
  ns = (int)Emax.size();
  L = mapper.getNormalisedSystemLength();
  Array1dr xj( mapper.getNormalisedGridPoints() );

  v.resize(ns);
  Emin.resize(ns);
  q.resize(ns);
  n0.resize(ns);
  u.resize(ns);
  for (int s=0; s<ns; ++s) {
    // Eq. 31
    v[s] = m[s]*4.0*m_pi/L;
    Emin[s] = 0.0;
    q[s] = 1.0;
    // Calculation of n0
    // Eq. 27
    Array1dr w(Emax[s]/sqrt(2.0*(1.0-pow2(v[s])))*xj);
    // Eq. 26
    Array1dr F(Emax[s]/cosh(w));
    // Eq. 22
    Array1dr n(-pow2(F)/(1.0-pow2(v[s])));
    n0[s] = -integrate1d(xj, n)/L;
    // Eq. 29
    u[s] = 0.5*v[s]+2.0*n0[s]/v[s]-
           (pow2(Emax[s])+pow2(Emin[s]))/(v[s]*(1.0-pow2(v[s])));
  }
}

void Solitons::initEndn(Mapper &mapper, FT1d &fft,
                        Array1dc &Ej, Array1dc &nj, Array1dc &dnj,
                        Array1dc &Ek, Array1dc &nk, Array1dc &dnk)
{

  InitMode::initEndn(mapper, fft, Ej, nj, dnj, Ek, nk, dnk);

  Array1dr xj( mapper.getNormalisedGridPoints() );

  for (int s=0; s<ns; ++s) {
    // Eq. 24
    Real Phi = v[s]/2.0;
    // Eq. 27
    Real t = 0.0;
    Array1dr x(xj+x0[s]-v[s]*t);
#ifdef BZ_HAVE_SYSTEM_V_MATH
    x = fmod(x+L/2.0, L)-L/2.0;
#else
    x -= L*floor((x+L/2.0)/L); // Calculate x modulo [-L/2, L/2]
#endif
    Array1dr w(Emax[s]/sqrt(2.0*(1.0-pow2(v[s])))*x);
    // Eq. 26
    Array1dr F(Emax[s]/cosh(w));
    x = xj+x0[s]-u[s]*t;
    Array1dc z(I*Phi*x);
    Array1dc E(-I*F*complexExp(z));
    // Eq. 22
    Array1dc n(zip(-pow2(F)/(1.0-pow2(v[s]))+n0[s], 0.0, Complex()));
    Array1dr dF(-v[s]*pow2(Emax[s])/sqrt(2.0*(1.0-pow2(v[s])))*
                sinh(w)/pow2(cosh(w)));
    Array1dc dn(zip(2.0/(1.0-pow2(v[s]))*F*dF, 0.0, Complex()));
    Ej += E;
    nj += n;
    dnj += dn;
  }

  fft.direct(Ej, Ek);
  fft.direct(nj, nk);
  fft.direct(dnj, dnk);

#if 1
  saveMatlab("solitons.m", ios::out|ios::trunc, "xj", xj);
  saveMatlab("solitons.m", ios::out|ios::app, "Ej", Ej);
  saveMatlab("solitons.m", ios::out|ios::app, "nj", nj);
  saveMatlab("solitons.m", ios::out|ios::app, "dnj", dnj);
  Array1dr k( mapper.getNormalisedWaveNumbers() );
  saveMatlab("solitons.m", ios::out|ios::app, "k", k);
  saveMatlab("solitons.m", ios::out|ios::app, "Ek", Ek);
  saveMatlab("solitons.m", ios::out|ios::app, "nk", nk);
  saveMatlab("solitons.m", ios::out|ios::app, "dnk", dnk);
#endif
}

void Solitons::printOn(ostream& os) const
{
  InitMode::printOn(os);
  os
      << '\n'
      << "Emax = " << Emax << '\n'
      << "Emin = " << Emin << '\n'
      << "m    = " << m << '\n'
      << "v    = " << v << '\n'
      << "u    = " << u << '\n'
      << "n0   = " << n0 << '\n'
      << "x0   = " << x0 << '\n'
      << "q    = " << q;
}

void Solitons::initParsing(int nargs, char *args[])
{
  const string id(PACKAGE);
  const string version(
    VERSION "\n"
    "$Id: solitons.cpp,v 1.48 2011/03/26 07:47:28 patrick Exp $"
    "\n");
  const string copyright(ZAK_COPYRIGHT);

  registerClass("Solitons");
  registerPackage(id, version, copyright);

  using parser::types::intVect;
  using parser::types::realVect;

  insertOption(_Emax, "Emax", realVect, "Emax of solitons", Any());
  insertOption(_m, "m", intVect, "Modulo m of solitons", Any());
  insertOption(_x0, "x0", realVect, "Offset x0 of solitons", Any());
}

void Solitons::parseParam()
{
  parseOption(_Emax, Emax);
  parseOption(_m, m);
  parseOption(_x0, x0);
}

void Solitons::checkParam() const
{
  if (Emax.size() != m.size() ||  Emax.size() != x0.size()) {
    ostringstream os;
    os << "Emax, m, x0 have different size\n"
       << " size(Emax) = " << Emax.size() << '\n'
       << " size(m)    = " << m.size() << '\n'
       << " size(x0)   = " << x0.size();
    throw ClassException("Solitons", os.str());
  }
}

