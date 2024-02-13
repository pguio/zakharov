/**************************************************************************
 *
 * $Id: egauss.h,v 1.26 2011/03/26 07:47:28 patrick Exp $
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
 **************************************************************************/

#ifndef EGAUSS_H
#define EGAUSS_H

#include <default-mode.h>

class EGauss: public InitMode {
public:

  typedef std::ostream ostream;

  friend ostream& operator<<(ostream& os, const EGauss &egauss);

  EGauss(int nargs, char *args[]);
  ~EGauss();

  virtual void initialise(Mapper &mapper);
  void initEndn(Mapper &mapper, FT1d &fft,
                Array1dc &Ej, Array1dc &nj, Array1dc &dnj,
                Array1dc &Ek, Array1dc &nk, Array1dc &dnk);

protected:

  enum parser_enum {
    _E0=InitMode::next, _Erms, _nrms, _k0, _W, _Delta
  };

  virtual void printOn(ostream& os) const;

private:
  Real E0, Erms, nrms, k0, W, Delta;

  void initParsing(int nargs, char *args[]);
  void parseParam();
};

#endif // EGAUSS_H

