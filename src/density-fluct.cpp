/**************************************************************************
 *
 * $Id: density-fluct.cpp,v 1.12 2020/04/22 13:34:53 patrick Exp $
 *
 * Copyright (c) 2010-2011 Patrick Guio <patrick.guio@gmail.com>
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
#include <density-fluct.h>

using std::ios;
using std::ostream;

namespace {
  RegisterInFactory<InitMode, DensityFluct> registerMe("densfluct");
}

namespace factory {
  void dummyDensityFluct()
  {}
}

ostream& operator<<(ostream& os, const DensityFluct &fluct)
{
  fluct.printOn(os);
  return os;
}

DensityFluct::DensityFluct(int nargs, char *args[]) :
  InitMode(nargs, args, "DensityFluct"), dn(0), l0(0), lm(0)
{
  initParsing(nargs, args);
  parseParam();
  checkParam();
}

DensityFluct::~DensityFluct()
{}

void DensityFluct::initialise(Mapper &mapper)
{
  InitMode::initialise(mapper);

#if 0
  source_nk(0) = 0.0;
  source_Ek(0) = 0.0;
#else
  source_nk = mapper.getnkFluctuationsLevel()*N;
  source_Ek = mapper.getEkFluctuationsLevel()*N;
#endif

  Array1dr k( mapper.getNormalisedWaveNumbers() );

  Range isaved(Range(0,adjfn(0,N-1,8),8));
  saveMatlab("tfluct.m", ios::out|ios::trunc, "k", k(isaved));
  saveMatlab("tfluct.m", ios::out|ios::app, "snk", source_nk(isaved));
  saveMatlab("tfluct.m", ios::out|ios::app, "sEk", source_Ek(isaved));
}

void DensityFluct::initEndn(Mapper &mapper, FT1d &fft,
                            Array1dc &Ej, Array1dc &nj, Array1dc &dnj,
                            Array1dc &Ek, Array1dc &nk, Array1dc &dnk)
{
  using std::cout;
  using std::endl;
  using blitz::pow2;
  using blitz::real;
  using blitz::max;
  using blitz::abs;

  // default init mode with thermal fluctuations
  InitMode::initEndn(mapper, fft, Ej, nj, dnj, Ek, nk, dnk);

  // physical parameters
  Array1dr k( fourier::ifftshift(mapper.getPhysicalWaveNumbers() ) );

  Real k0 = m_2pi / l0;
  Real km = m_2pi / lm;
  Real dk = km-k0;

  Real ne = mapper.getDensity();
  dnOverne = dn/ne;

  Array1dc densfluctk(N), phase(N);

  densfluctk = Complex(0.0, 0.0);
  phase = Complex(0.0, 0.0);

  for (int i=1; i<M; ++i) {
    phase(i)        = Complex(0.0, m_2pi*rngn.random());
    phase(N-i)      = -phase(i);
    densfluctk(i)   = exp(-pow2((k(i)-k0)/dk));
    densfluctk(N-i) = exp(-pow2((k(N-i)+k0)/dk));
  }
  densfluctk *= complexExp(phase);

  Array1dc densfluctj(N);
  fft.inverse(densfluctk, densfluctj);
  // to physical density
  densfluctj = mapper.toPhysicalDensity(densfluctj);

  // scale dn/ne
  Real mxdensfluct = max(abs(real(densfluctj)));
  densfluctj = dn*densfluctj/mxdensfluct;
  mxdensfluct = max(abs(real(densfluctj)));
  densfluctj = blitz::zip(real(densfluctj), 0.0, Complex());
  cout << "mxdensfluct = " << mxdensfluct << " [m-3], "
       << "mxdensfluct/ne = " << mxdensfluct/ne*100 << " [%]" << endl;

  // back to normalised units
  densfluctj = mapper.toNormalisedDensity(densfluctj);

  // Add density fluctuation to thermal fluctuations
  nj += densfluctj;

  // Ditto k-spectrum
  fft.direct(densfluctj, densfluctk);
  nk += densfluctk;

  source_nk = densfluctk;

  saveMatlab("densfluct.m", ios::out|ios::trunc, "dn", dn);
  saveMatlab("densfluct.m", ios::out|ios::app, "dnOverne", dnOverne);
  saveMatlab("densfluct.m", ios::out|ios::app, "l0", l0);
  saveMatlab("densfluct.m", ios::out|ios::app, "lm", lm);
#if 0
  Range isaved(Range(0,adjfn(0,N-1,8),8));
  saveMatlab("densfluct.m", ios::out|ios::trunc, "k", k(isaved));
  saveMatlab("densfluct.m", ios::out|ios::app, "source_nk", source_nk(isaved));
#endif
}

void DensityFluct::printOn(ostream& os) const
{
  InitMode::printOn(os);
  os
      << '\n'
      << "dn    = " << dn << " [m-3]\n"
      << "dn/ne = " << dnOverne*100 << " [%]\n"
      << "l0    = " << l0 << " [m]\n"
      << "lm    = " << lm << " [m]";
}

void DensityFluct::initParsing(int nargs, char *args[])
{
  const string id(PACKAGE);
  const string version(
    VERSION "\n"
    "$Id: density-fluct.cpp,v 1.12 2020/04/22 13:34:53 patrick Exp $"
    "\n");
  const string copyright(ZAK_COPYRIGHT);

  registerClass("DensityFluct");
  registerPackage(id, version, copyright);

  using parser::types::real;

  insertOption(_dn, "dn", real, "Parameter delta n [m-3]", Any(dn));
  insertOption(_l0, "l0", real, "Parameter lambda_0 [m]", Any(l0));
  insertOption(_lm, "lm", real, "Parameter lambda_min [m]", Any(lm));
}

void DensityFluct::parseParam()
{

  parseOption(_dn, dn);
  parseOption(_l0, l0);
  parseOption(_lm, lm);

}

void DensityFluct::checkParam() const
{
  if (l0 != 0.0 && l0 <= lm) {
    throw ClassException("DensityFluct","lm must be strictly smaller than l0");
  }
}

