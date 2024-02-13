/**************************************************************************
 *
 * $Id: default-mode.h,v 1.26 2011/03/26 07:47:28 patrick Exp $
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

#ifndef DEFAULT_MODE_H
#define DEFAULT_MODE_H

#include <fstream>
#include <iostream>

#include <zakharov-defs.h>
#include <mapper.h>

class InitMode : public parser::Parser {
public:

  typedef std::ostream ostream;
  typedef std::string string;

  typedef blitz::Range Range;

  typedef fourier::ODFT1D<fourier::complex,fourier::complex> FT1d;

  friend ostream& operator<<(ostream& os, const InitMode &mode);

  InitMode(int nargs, char *args[], const string modename="DefaultMode");
  virtual ~InitMode();

  virtual void initialise(Mapper &mapper);
  virtual void initEndn(Mapper &mapper, FT1d &fft,
                        Array1dc &Ej, Array1dc &nj, Array1dc &dnj,
                        Array1dc &Ek, Array1dc &nk, Array1dc &dnk);

  void getSources(Array1dc &delta_nk, Array1dc &delta_Ek);

protected:

  typedef parser::Parser Parser;

  enum parser_enum {
    _mseed=1, next
  };

  const string ModeName;
  Vector2di mseed; // mseed(0) -> rngE, mseed(1) -> rngn
  Array1dc source_nk;
  Array1dc source_Ek;
  int N, M;

  typedef ranlib::MersenneTwister MersenneTwister;
  typedef ranlib::sharedState sharedState; // shares an IRNG with other RNGs
  typedef ranlib::independentState independentState; // contains its own IRNG
  ranlib::UniformOpenClosed<Real, MersenneTwister, independentState> rngE;
  ranlib::UniformOpenClosed<Real, MersenneTwister, independentState> rngn;

  virtual void printOn(ostream& os) const;

private:

  void initParsing(int nargs, char *args[]);
  void parseParam();
};

#endif // DEFAULT_MODE_H
