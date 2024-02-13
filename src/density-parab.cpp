/**************************************************************************
 *
 * $Id: density-parab.cpp,v 1.5 2020/04/22 13:34:53 patrick Exp $
 *
 * Copyright (c) 2011 Patrick Guio <patrick.guio@gmail.com>
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
#include <density-parab.h>

using std::ios;
using std::ostream;

namespace {
  RegisterInFactory<InitMode, DensityParab> registerMe("densparab");
}

namespace factory {
  void dummyDensityParab()
  {}
}

ostream& operator<<(ostream& os, const DensityParab &fluct)
{
  fluct.printOn(os);
  return os;
}

DensityParab::DensityParab(int nargs, char *args[]) :
  InitMode(nargs, args, "DensityParab"), dn(0), dl(0)
{
  initParsing(nargs, args);
  parseParam();
}

DensityParab::~DensityParab()
{}

void DensityParab::initialise(Mapper &mapper)
{
  InitMode::initialise(mapper);

  source_nk = mapper.getnkFluctuationsLevel()*N;
  source_Ek = mapper.getEkFluctuationsLevel()*N;
}

void DensityParab::initEndn(Mapper &mapper, FT1d &fft,
                            Array1dc &Ej, Array1dc &nj, Array1dc &dnj,
                            Array1dc &Ek, Array1dc &nk, Array1dc &dnk)
{
  using std::cout;
  using std::endl;
  using blitz::pow2;
  using blitz::real;
  using blitz::max;
  using blitz::abs;

  InitMode::initEndn(mapper, fft, Ej, nj, dnj, Ek, nk, dnk);

  // physical parameters
  Array1dr x( mapper.getPhysicalGridPoints() );

  Real ne = mapper.getDensity();
  dnOverne = fabs(dn)/ne;

  // from Ergun et al, 2008 delta n = n_0 K0^2 x^2/2
  K0 = sqrt(2.0*dnOverne)/(dl/2.0);

  Array1dc densparabj(N), densparabk(N);

  if ( dn > 0.0)
    densparabj = where(abs(x)<dl/2,
                       ne * pow2(K0) * pow2(x) / 2.0,
                       ne * pow2(K0) * pow2(dl/2)/ 2.0);
  else
    densparabj = where(abs(x)<dl/2,
                       - ne * pow2(K0) * pow2(x) / 2.0,
                       - ne * pow2(K0) * pow2(dl/2)/ 2.0);

  // convert to normalised parameters
  densparabj = mapper.toNormalisedDensity(densparabj);

  fft.direct(densparabj, densparabk);
  // dn(k=0) = 0
  densparabk(0) = 0.0;
  // zero-padding
  densparabk(Range(M+1,(N-1)-M)) = 0.0;
  fft.inverse(densparabk, densparabj);

  nj += densparabj;
  nk += densparabk;

  // source term spectra identical to initial spectra
  source_nk = densparabk;

  saveMatlab("densparab.m", ios::out|ios::trunc, "dn", dn);
  saveMatlab("densparab.m", ios::out|ios::app, "dl", dl);
  saveMatlab("densparab.m", ios::out|ios::app, "K0", K0);
}

void DensityParab::printOn(ostream& os) const
{
  InitMode::printOn(os);
  os
      << '\n'
      << "dn    = " << dn << " [m-3]\n"
      << "dl    = " << dl << " [m]\n"
      << "dn/ne = " << dnOverne*100 << " [%]\n"
      << "K0    = " << K0 << " [m-2]";
}

void DensityParab::initParsing(int nargs, char *args[])
{
  const string id(PACKAGE);
  const string version(
    VERSION "\n"
    "$Id: density-parab.cpp,v 1.5 2020/04/22 13:34:53 patrick Exp $"
    "\n");
  const string copyright(ZAK_COPYRIGHT);

  registerClass("DensityParab");
  registerPackage(id, version, copyright);

  using parser::types::real;

  insertOption(_dn, "dn", real, "Parameter delta n [m-3]", Any(dn));
  insertOption(_dl, "dl", real, "Parameter delta l [m]", Any(dl));
}

void DensityParab::parseParam()
{

  parseOption(_dn, dn);
  parseOption(_dl, dl);

}
