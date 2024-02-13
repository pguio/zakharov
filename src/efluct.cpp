/**************************************************************************
 *
 * $Id: efluct.cpp,v 1.31 2011/03/26 07:47:28 patrick Exp $
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
#include <efluct.h>

using std::ostream;

namespace {
  RegisterInFactory<InitMode, EFluct> registerMe("efluct");
}

namespace factory {
  void dummyEfluct()
  {}
}


ostream& operator<<(ostream& os, const EFluct &e)
{
  e.printOn(os);
  return os;
}

EFluct::EFluct(int nargs, char *args[]) : InitMode(nargs, args, "EFluct"),
  E0(0.1), Erms(1.0e-4)
{
  initParsing(nargs, args);
  parseParam();
}

EFluct::~EFluct()
{}

void EFluct::initialise(Mapper &mapper)
{
  InitMode::initialise(mapper);

  source_Ek(0) = E0;
  source_Ek(Range(1,N-1)) = Erms;

  source_Ek *= (Real)N;
}

void EFluct::printOn(ostream& os) const
{
  InitMode::printOn(os);
  os
      << '\n'
      << "E0   = " << E0 << '\n'
      << "Erms = " << Erms;
}

void EFluct::initParsing(int nargs, char *args[])
{
  const string id(PACKAGE);
  const string version(
    VERSION "\n"
    "$Id: efluct.cpp,v 1.31 2011/03/26 07:47:28 patrick Exp $"
    "\n");
  const string copyright(ZAK_COPYRIGHT);

  registerClass("EFluct");
  registerPackage(id, version, copyright);

  using parser::types::real;

  insertOption(_E0  , "E0"  , real, "Parameter E0 " , Any(E0));
  insertOption(_Erms, "Erms", real, "Parameter Erms", Any(Erms));
}

void EFluct::parseParam()
{
  parseOption(_E0  , E0);
  parseOption(_Erms, Erms);
}
