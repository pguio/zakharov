/**************************************************************************
 *
 * $Id: solitons.h,v 1.27 2011/03/26 07:47:28 patrick Exp $
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

#ifndef SOLITONS_H
#define SOLITONS_H

#include <default-mode.h>

class Solitons : public InitMode {
public:

  typedef std::ostream ostream;
  typedef std::ostringstream ostringstream;

  friend ostream& operator<<(ostream& os, const Solitons &solitons);

  Solitons(int nargs, char *args[]);
  virtual ~Solitons();

  virtual void initialise(Mapper &mapper);
  void initEndn(Mapper &mapper, FT1d &fft,
                Array1dc &Ej, Array1dc &nj, Array1dc &dnj,
                Array1dc &Ek, Array1dc &nk, Array1dc &dnk);

protected:

  enum parser_enum {
    _Emax=InitMode::next, _m, _x0
  };

  virtual void printOn(ostream& os) const;

private:

  typedef std::vector<Real> vectorr;
  typedef std::vector<int> vectori;

  vectorr Emax;
  vectori m;
  vectorr x0;
  vectorr Emin;
  vectorr v;
  vectorr Phi;
  vectorr n0;
  vectorr u;
  vectorr q;
  Real L;
  int ns;

  void initParsing(int nargs, char *args[]);
  void parseParam();
  void checkParam() const;
};

#endif // SOLITONS_H

