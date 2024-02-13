/**************************************************************************
 *
 * $Id: thermal-fluct.cpp,v 1.40 2011/03/26 07:47:28 patrick Exp $
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
#include <thermal-fluct.h>

using std::ios;
using std::ostream;

namespace {
  RegisterInFactory<InitMode, ThermalFluct> registerMe("tfluct");
}

namespace factory {
  void dummyThermalFluct()
  {}
}

ostream& operator<<(ostream& os, const ThermalFluct &fluct)
{
  fluct.printOn(os);
  return os;
}

ThermalFluct::ThermalFluct(int nargs, char *args[]) :
  InitMode(nargs, args, "ThermalFluct")
{
  initParsing(nargs, args);
  parseParam();
}

ThermalFluct::~ThermalFluct()
{}

void ThermalFluct::initialise(Mapper &mapper)
{
  InitMode::initialise(mapper);

  source_nk = mapper.getnkFluctuationsLevel()*N;
  source_Ek = mapper.getEkFluctuationsLevel()*N;
#if 0
  source_nk(0) = 0.0;
  source_Ek(0) = 0.0;
#endif
  Array1dr k( mapper.getNormalisedWaveNumbers() );

  Range isaved(Range(0,adjfn(0,N-1,8),8));
  saveMatlab("tfluct.m", ios::out|ios::trunc, "k", k(isaved));
  saveMatlab("tfluct.m", ios::out|ios::app, "snk", source_nk(isaved));
  saveMatlab("tfluct.m", ios::out|ios::app, "sEk", source_Ek(isaved));
}

void ThermalFluct::printOn(ostream& os) const
{
  InitMode::printOn(os);
}

void ThermalFluct::initParsing(int nargs, char *args[])
{
  const string id(PACKAGE);
  const string version(
    VERSION "\n"
    "$Id: thermal-fluct.cpp,v 1.40 2011/03/26 07:47:28 patrick Exp $"
    "\n");
  const string copyright(ZAK_COPYRIGHT);

  registerClass("ThermalFluct");
  registerPackage(id, version, copyright);
}

void ThermalFluct::parseParam()
{}

