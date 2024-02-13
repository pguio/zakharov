/**************************************************************************
 *
 * $Id: density-fluct.h,v 1.4 2011/03/26 07:47:28 patrick Exp $
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

#ifndef DENSITY_FLUCT_H
#define DENSITY_FLUCT_H

#include <default-mode.h>

class DensityFluct: public InitMode {
public:

  friend ostream& operator<<(ostream& os, const DensityFluct &fluct);

  DensityFluct(int nargs, char *args[]);
  ~DensityFluct();

  void initialise(Mapper &mapper);
  void initEndn(Mapper &mapper, FT1d &fft,
                Array1dc &Ej, Array1dc &nj, Array1dc &dnj,
                Array1dc &Ek, Array1dc &nk, Array1dc &dnk);

protected:

  enum parser_enum {
    _dn = InitMode::next, _l0, _lm
  };

  void printOn(ostream& os) const;

private:

  Real dn, l0, lm;
  Real dnOverne;

  void initParsing(int nargs, char *args[]);
  void parseParam();
  void checkParam() const;

};

#endif // DENSITY_FLUCT_H

