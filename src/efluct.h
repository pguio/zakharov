/**************************************************************************
 *
 * $Id: efluct.h,v 1.23 2011/03/26 07:47:28 patrick Exp $
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

#ifndef EFLUCT_H
#define EFLUCT_H

#include <default-mode.h>

class EFluct: public InitMode {
public:

  typedef std::ostream ostream;

  typedef blitz::Range Range;

  friend ostream& operator<<(ostream& os, const EFluct &e);

  EFluct(int nargs, char *args[]);
  ~EFluct();

  virtual void initialise(Mapper &mapper);

protected:

  enum parser_enum {
    _E0=InitMode::next, _Erms
  };

  virtual void printOn(ostream& os) const;

private:

  Real E0, Erms;

  void initParsing(int nargs, char *args[]);
  void parseParam();

};

#endif // EFLUCT_H

