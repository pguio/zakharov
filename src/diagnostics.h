/**************************************************************************
 *
 * $Id: diagnostics.h,v 1.66 2011/03/26 07:47:28 patrick Exp $
 *
 * Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
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

#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include <fstream>
#include <iostream>

#include <zakharov-defs.h>
#include <mapper.h>
#include <spectra.h>
#include <hdf-interface.h>

class Diagnostics : public parser::Parser {
public:

  typedef std::string string;
  typedef std::ostream ostream;

  friend ostream& operator<<(ostream& os, const Diagnostics &d);

  Diagnostics(int nargs, char *args[]);
  ~Diagnostics();

  bool parseHelp() const;
  bool parseVersion() const;
  bool parseTemplate() const;

  void startInit();
  void initnE(Mapper &mapper);
  void initNPH(Mapper &mapper);
  void initR(Mapper &mapper);
  void initW(Mapper &mapper);
  void initDiffusion(Mapper &mapper);
  void initSk2(Mapper &mapper, Spectra &spectra);
  void initSkf2(Mapper &mapper, Spectra &spectra);
  void endInit();

  bool saveStandardDiagnostics(int i) const {
    return iter.isWithin(i) && ((i-iter.start()) % iter.stride()) == 0;
  }

  void savenE(Mapper &mapper, int i,
              Array1dc &nj, Array1dc &Ej, Array1dc &nk, Array1dc &Ek);
  void saveNPH(Mapper &mapper, int i, Real N, Real P, Real H);
  void saveR(Mapper &mapper, int i, Real R);
  void saveW(Mapper &mapper, int i, Real W);
  void saveDiffusion(Mapper &mapper, int i);
  void saveSk2(Mapper &mapper, Spectra &spectra);
  void saveSkf2(Mapper &mapper, int i, bool atEnd, Spectra &spectra);

  string getDiagnosticFilename() const {
    return hdf.getFilename();
  }

protected:

  typedef parser::Parser Parser;

private:

  enum parser_enum { _nj=1, _Ej, _nk, _Ek, _NPH, _iter, _Fe, _nue,
                     _Sk2, _Skf2, _Seef
                   };
  enum dim_enum { time, x, y, z};

  RangeSpec<int> iter;
  bool njDiag, EjDiag, nkDiag, EkDiag, NPHDiag,
       FeDiag, nueDiag, Sk2Diag, Skf2Diag, SeefDiag;

  hdf::Hdf hdf;
  hdf::SDvar V;
  hdf::VecInt start, edge;

  void initParsing(int nargs, char *args[]);
  void paramParsing();

  void initVector(const string &name, Real dt);
  void init1dField(const string &name, Array1dr &axis, const string &axis_name);
  void init2dField(const string &name, Real dt,
                   const string &axis_name, Array1dr &axis);
  void init2dField(const string &name,
                   Array1dr &axis1, const string &axis_name1,
                   Array1dr &axis2, const string &axis_name2);
  void init3dField(const string &name,
                   Array1dr &axis1, const string &axis_name1,
                   Array1dr &axis2, const string &axis_name2,
                   Array1dr &axis3, const string &axis_name3);
};


#endif // DIAGNOSTICS_H
