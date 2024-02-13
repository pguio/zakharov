/**************************************************************************
 *
 * $Id: diagnostics.cpp,v 1.86 2014/10/16 13:24:28 patrick Exp $
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

#include <diagnostics.h>

using blitz::imag;
using blitz::real;

using parser::header;
using parser::yesno;

using std::cout;
using std::endl;
using std::ostream;

using std::string;

const string re_nj_hdf_name("re-nj");
const string im_nj_hdf_name("im-nj");
const string re_Ej_hdf_name("re-Ej");
const string im_Ej_hdf_name("im-Ej");
const string xj_hdf_name("xj");

const string re_nk_hdf_name("re-nk");
const string im_nk_hdf_name("im-nk");
const string re_Ek_hdf_name("re-Ek");
const string im_Ek_hdf_name("im-Ek");
const string k_hdf_name("k");

const string N_hdf_name("N");
const string P_hdf_name("P");
const string H_hdf_name("H");
const string R_hdf_name("R");
const string W_hdf_name("W");

const string Fe_hdf_name("Fe");
const string Fe_v_hdf_name("v-Fe");
const string nue_hdf_name("nue");

const string nk2_hdf_name("nk2");
const string Ek2_hdf_name("Ek2");
const string Sk2_k_hdf_name("k-Sk2");
const string Sk2_t_hdf_name("t-Sk2");

const string nkf2_hdf_name("nkf2");
const string nkf2_f_hdf_name("f-nkf2");
const string nkf2_k_hdf_name("k-nkf2");
const string nkf2_t_hdf_name("t-nkf2");

const string Ekf2_hdf_name("Ekf2");
const string Ekf2_f_hdf_name("f-Ekf2");
const string Ekf2_k_hdf_name("k-Ekf2");
const string Ekf2_t_hdf_name("t-Ekf2");

const string ukf2_hdf_name("ukf2");
const string ukf2_f_hdf_name("f-ukf2");
const string ukf2_k_hdf_name("k-ukf2");
const string ukf2_t_hdf_name("t-ukf2");

const string dkf2_hdf_name("dkf2");
const string dkf2_f_hdf_name("f-dkf2");
const string dkf2_k_hdf_name("k-dkf2");
const string dkf2_t_hdf_name("t-dkf2");

const string seef_hdf_name("seef");
const string seef_f_hdf_name("f-seef");
const string seef_t_hdf_name("t-seef");


const string time_hdf_name("time");

ostream& operator<<(ostream& os, const Diagnostics &d)
{
  return os
         << d.hdf << '\n'
         << header("Diagnostics setup")
         << "iter      = " << d.iter << '\n'
         << "nj        = " << yesno(d.njDiag) << '\n'
         << "Ej        = " << yesno(d.EjDiag) << '\n'
         << "nk        = " << yesno(d.nkDiag) << '\n'
         << "Ek        = " << yesno(d.EkDiag) << '\n'
         << "NPH       = " << yesno(d.NPHDiag) << '\n'
         << "Fe        = " << yesno(d.FeDiag) << '\n'
         << "nue       = " << yesno(d.nueDiag) << '\n'
         << "Sk2       = " << yesno(d.Sk2Diag) << '\n'
         << "Skf2      = " << yesno(d.Skf2Diag);
}


Diagnostics::Diagnostics(int nargs, char *args[]) : Parser(nargs, args),
  iter(DEFAULT_SaveIter),
  njDiag(false), EjDiag(false), nkDiag(false), EkDiag(false), NPHDiag(true),
  FeDiag(false), nueDiag(false), Sk2Diag(false), Skf2Diag(false),
  hdf(nargs, args)
{
  initParsing(nargs, args);
  paramParsing();
}

Diagnostics::~Diagnostics()
{}


bool Diagnostics::parseHelp() const
{
  if (Parser::parseHelp()) {
    hdf.parseHelp();
    return true;
  }
  return false;
}

bool Diagnostics::parseVersion() const
{
  if (Parser::parseVersion() ) {
    hdf.parseVersion();
    return true;
  }
  return false;
}

bool Diagnostics::parseTemplate() const
{
  if (Parser::parseTemplate() ) {
    hdf.parseTemplate();
    return true;
  }
  return false;
}

void Diagnostics::initParsing(int nargs, char *args[])
{
  const string id(PACKAGE);
  const string version(
    VERSION "\n"
    "$Id: diagnostics.cpp,v 1.86 2014/10/16 13:24:28 patrick Exp $"
    "\n");
  const string copyright(ZAK_COPYRIGHT);

  registerClass("Diagnostics");
  registerPackage(id, version, copyright);

  using parser::types::boolean;
  using parser::types::intVect;

  insertOption(_nj  , "nj"  , boolean, "Diagnostic nj"                              , Any(njDiag));
  insertOption(_Ej  , "Ej"  , boolean, "Diagnostic Ej"                              , Any(EjDiag));
  insertOption(_nk  , "nk"  , boolean, "Diagnostic nk"                              , Any(nkDiag));
  insertOption(_Ek  , "Ek"  , boolean, "Diagnostic Ek"                              , Any(EkDiag));
  insertOption(_NPH , "NPH" , boolean, "Diagnostic NPH"                             , Any(NPHDiag));
  insertOption(_iter, "iter", intVect, "Diagnostic save (start,stride,end)"         , Any(iter));
  insertOption(_Fe  , "Fe"  , boolean, "Diagnostic Fe"                              , Any(FeDiag));
  insertOption(_nue , "nue" , boolean, "Diagnostic nue"                             , Any(nueDiag));
  insertOption(_Sk2 , "Sk2" , boolean, "Diagnostic <|E(k)|^2> and  <|n(k)|^2>"      , Any(Sk2Diag));
  insertOption(_Skf2, "Skf2", boolean, "Diagnostic k^2<|E(k,f)|^2> and <|n(k,f)|^2>", Any(Skf2Diag));
  insertOption(_Seef, "Seef", boolean, "Diagnostic SEE spectrum <|n(k=0,f)E(k=0,f)|> ", Any(SeefDiag));
}

void Diagnostics::paramParsing()
{
  parseOption(_nj  , njDiag);
  parseOption(_Ej  , EjDiag);
  parseOption(_nk  , nkDiag);
  parseOption(_Ek  , EkDiag);
  parseOption(_NPH , NPHDiag);
  parseOption(_iter, iter.data);
  parseOption(_Fe  , FeDiag);
  parseOption(_nue , nueDiag);
  parseOption(_Sk2 , Sk2Diag);
  parseOption(_Skf2, Skf2Diag);
  parseOption(_Seef, SeefDiag);
}


void Diagnostics::startInit()
{
  if (njDiag || EjDiag || nkDiag || EkDiag || NPHDiag ||
      FeDiag || nueDiag || Sk2Diag || Skf2Diag || SeefDiag)
    hdf.createSD(
      PACKAGE " " VERSION
      " $Id: diagnostics.cpp,v 1.86 2014/10/16 13:24:28 patrick Exp $");
}

void Diagnostics::initnE(Mapper &mapper)
{
  // First dimension: time
  Real dt = mapper.getPhysicalTimeStep();

  // Second dimension: physical space
  Array1dr xj( mapper.getPhysicalGridPoints() );
  if (njDiag) {
    init2dField(re_nj_hdf_name, dt, xj_hdf_name, xj);
    hdf.initSDvar(V);
    init2dField(im_nj_hdf_name, dt, xj_hdf_name, xj);
    hdf.initSDvar(V);
  }
  if (EjDiag) {
    init2dField(re_Ej_hdf_name, dt, xj_hdf_name, xj);
    hdf.initSDvar(V);
    init2dField(im_Ej_hdf_name, dt, xj_hdf_name, xj);
    hdf.initSDvar(V);
  }

  // Second dimension: Fourier wavenumber
  Array1dr k ( mapper.getPhysicalWaveNumbers() );
#if 1
  // Remove zero-padded wavenumbers
  int N = mapper.getGridSize();
  int M = mapper.getZeroPaddingSize();
  Array1dr tmp(1.0*k);
  k.resize(2*M+1);
  k = tmp(blitz::Range(N/2-M,N/2+M));
#endif
  if (nkDiag) {
    init2dField(re_nk_hdf_name, dt, k_hdf_name, k);
    hdf.initSDvar(V);
    init2dField(im_nk_hdf_name, dt, k_hdf_name, k);
    hdf.initSDvar(V);
  }
  if (EkDiag) {
    init2dField(re_Ek_hdf_name, dt, k_hdf_name, k);
    hdf.initSDvar(V);
    init2dField(im_Ek_hdf_name, dt, k_hdf_name, k);
    hdf.initSDvar(V);
  }
}

void Diagnostics::initNPH(Mapper &mapper)
{
  Real dt = mapper.getPhysicalTimeStep();
  if (NPHDiag) {
    initVector(N_hdf_name, dt);
    hdf.initSDvar(V);
    initVector(P_hdf_name, dt);
    hdf.initSDvar(V);
    initVector(H_hdf_name, dt);
    hdf.initSDvar(V);
  }
}

void Diagnostics::initR(Mapper &mapper)
{
  Real dt = mapper.getPhysicalTimeStep();
  if (NPHDiag) {
    initVector(R_hdf_name, dt);
    hdf.initSDvar(V);
  }
}

void Diagnostics::initW(Mapper &mapper)
{
  Real dt = mapper.getPhysicalTimeStep();
  if (NPHDiag) {
    initVector(W_hdf_name, dt);
    hdf.initSDvar(V);
  }
}

void Diagnostics::initDiffusion(Mapper &mapper)
{
  Real dt = mapper.getPhysicalTimeStep();
  // Second dimension: Fourier wavenumber
  Array1dr k ( mapper.getPhysicalWaveNumbers() );
#if 1
  // Remove zero-padded wavenumbers
  int N = mapper.getGridSize();
  int M = mapper.getZeroPaddingSize();
  Array1dr tmp(1.0*k);
  k.resize(2*M+1);
  k = tmp(blitz::Range(N/2-M,N/2+M));
#endif
  if (FeDiag) {
    init2dField(Fe_hdf_name, dt, k_hdf_name, k);
    hdf.initSDvar(V);
  }
  if (nueDiag) {
    init2dField(nue_hdf_name, dt, k_hdf_name, k);
    hdf.initSDvar(V);
  }
}

void Diagnostics::initSk2(Mapper &mapper, Spectra &spectra)
{
  if (spectra.isKspecDiagnosed() && Sk2Diag) {
    Array1dr t(spectra.getKspecT());
    Array1dr k(spectra.getKspecK());
    init2dField(Ek2_hdf_name, t, Sk2_t_hdf_name, k, Sk2_k_hdf_name);
    hdf.initSDvar(V);
    init2dField(nk2_hdf_name, t, Sk2_t_hdf_name, k, Sk2_k_hdf_name);
    hdf.initSDvar(V);
  }
}

void Diagnostics::initSkf2(Mapper &mapper, Spectra &spectra)
{
  if (spectra.isIspecDiagnosed() && Skf2Diag) {
    Array1dr k(spectra.getIspecK());
    Array1dr f(spectra.getIspecF());
    if (spectra.isIspecIntegrated()) {
      init2dField(nkf2_hdf_name, k, nkf2_k_hdf_name, f, nkf2_f_hdf_name);
      hdf.initSDvar(V);
    } else {
      Array1dr t(spectra.getIspecT());
      init3dField(nkf2_hdf_name, t, nkf2_t_hdf_name, k, nkf2_k_hdf_name,
                  f, nkf2_f_hdf_name);
      hdf.initSDvar(V);
    }
  }
  if (spectra.isEspecDiagnosed() && Skf2Diag) {
    Array1dr k(spectra.getEspecK());
    Array1dr f(spectra.getEspecF());
    if (spectra.isEspecIntegrated()) {
      init2dField(Ekf2_hdf_name, k, Ekf2_k_hdf_name, f, Ekf2_f_hdf_name);
      hdf.initSDvar(V);
    } else {
      Array1dr t(spectra.getEspecT());
      init3dField(Ekf2_hdf_name, t, Ekf2_t_hdf_name, k, Ekf2_k_hdf_name,
                  f, Ekf2_f_hdf_name);
      hdf.initSDvar(V);
    }
  }
  if (spectra.isUspecDiagnosed() && Skf2Diag) {
    Array1dr k(spectra.getUspecK());
    Array1dr f(spectra.getUspecF());
    if (spectra.isUspecIntegrated()) {
      init2dField(ukf2_hdf_name, k, ukf2_k_hdf_name, f, ukf2_f_hdf_name);
      hdf.initSDvar(V);
    } else {
      Array1dr t(spectra.getUspecT());
      init3dField(ukf2_hdf_name, t, ukf2_t_hdf_name, k, ukf2_k_hdf_name,
                  f, ukf2_f_hdf_name);
      hdf.initSDvar(V);
    }
  }
  if (spectra.isDspecDiagnosed() && Skf2Diag) {
    Array1dr k(spectra.getDspecK());
    Array1dr f(spectra.getDspecF());
    if (spectra.isDspecIntegrated()) {
      init2dField(dkf2_hdf_name, k, dkf2_k_hdf_name, f, dkf2_f_hdf_name);
      hdf.initSDvar(V);
    } else {
      Array1dr t(spectra.getDspecT());
      init3dField(dkf2_hdf_name, t, dkf2_t_hdf_name, k, dkf2_k_hdf_name,
                  f, dkf2_f_hdf_name);
      hdf.initSDvar(V);
    }
  }
  if (spectra.isSEEspecDiagnosed() && SeefDiag) {
    Array1dr f(spectra.getSEEspecF());
    if (spectra.isSEEspecIntegrated()) {
      init1dField(seef_hdf_name, f, seef_f_hdf_name);
      hdf.initSDvar(V);
    } else {
      Array1dr t(spectra.getSEEspecT());
      init2dField(seef_hdf_name, t, seef_t_hdf_name, f, seef_f_hdf_name);
      hdf.initSDvar(V);
    }
  }
}

void Diagnostics::endInit()
{
#if 0
  if (njDiag || EjDiag || nkDiag || EkDiag || NPHDiag ||
      FeDiag || nueDiag || Sk2Diag || Skf2Diag)
    hdf.endInitSD();
#endif
}


void Diagnostics::initVector(const string &name, Real dt)
{
  V.resize(1);

  V.vartype =  hdf::type<double>::def;
  V.varname = name;

  V.dimsize[0] = iter.size();
  V.dimtype[0] = hdf::type<float>::def;
  V.dimname[0] = time_hdf_name;
  V.start[0] = dt*iter.start();
  V.stride[0] = dt*iter.stride();
  V.end[0] = dt*iter.end();
}

void Diagnostics::init1dField(const string &name,
                              Array1dr &axis, const string &axis_name)
{
  V.resize(1);

  V.varname = name;
  V.vartype =  hdf::type<double>::def;

  V.dimsize[0] = axis.rows();
  V.dimtype[0] = hdf::type<float>::def;
  V.dimname[0] = axis_name;
  V.start[0] = axis(0);
  V.stride[0] = (axis(V.dimsize[0]-1)-axis(0))/(V.dimsize[0]-1);
  V.end[0] = axis(V.dimsize[0]-1);
}

void Diagnostics::init2dField(const string &name, Real dt,
                              const string &axis_name, Array1dr &axis)
{
  V.resize(2);

  V.varname = name;
  V.vartype =  hdf::type<double>::def;

  V.dimsize[0] = iter.size();
  V.dimtype[0] = hdf::type<float>::def;
  V.dimname[0] = time_hdf_name;
  V.start[0] = dt*iter.start();
  V.stride[0] = dt*iter.stride();
  V.end[0] = dt*iter.end();

  V.dimsize[1] = axis.rows();
  V.dimtype[1] = hdf::type<float>::def;
  V.dimname[1] = axis_name;
  V.start[1] = axis(0);
  V.stride[1] = (axis(V.dimsize[1]-1)-axis(0))/(V.dimsize[1]-1);
  V.end[1] = axis(V.dimsize[1]-1);
}

void Diagnostics::init2dField(const string &name,
                              Array1dr &axis1, const string &axis_name1,
                              Array1dr &axis2, const string &axis_name2)
{
  V.resize(2);

  V.varname = name;
  V.vartype =  hdf::type<double>::def;

  V.dimsize[0] = axis1.rows();
  V.dimtype[0] = hdf::type<float>::def;
  V.dimname[0] = axis_name1;
  V.start[0] = axis1(0);
  V.stride[0] = (axis1(V.dimsize[0]-1)-axis1(0))/(V.dimsize[0]-1);
  V.end[0] = axis1(V.dimsize[0]-1);

  V.dimsize[1] = axis2.rows();
  V.dimtype[1] = hdf::type<float>::def;
  V.dimname[1] = axis_name2;
  V.start[1] = axis2(0);
  V.stride[1] = (axis2(V.dimsize[1]-1)-axis2(0))/(V.dimsize[1]-1);
  V.end[1] = axis2(V.dimsize[1]-1);
}

void Diagnostics::init3dField(const string &name,
                              Array1dr &axis1, const string &axis_name1,
                              Array1dr &axis2, const string &axis_name2,
                              Array1dr &axis3, const string &axis_name3)
{
  V.resize(3);

  V.varname = name;
  V.vartype =  hdf::type<double>::def;

  V.dimsize[0] = axis1.rows();
  V.dimtype[0] = hdf::type<float>::def;
  V.dimname[0] = axis_name1;
  V.start[0] = axis1(0);
  V.stride[0] = (axis1(V.dimsize[0]-1)-axis1(0))/(V.dimsize[0]-1);
  V.end[0] = axis1(V.dimsize[0]-1);

  V.dimsize[1] = axis2.rows();
  V.dimtype[1] = hdf::type<float>::def;
  V.dimname[1] = axis_name2;
  V.start[1] = axis2(0);
  V.stride[1] = (axis2(V.dimsize[1]-1)-axis2(0))/(V.dimsize[1]-1);
  V.end[1] = axis2(V.dimsize[1]-1);

  V.dimsize[2] = axis3.rows();
  V.dimtype[2] = hdf::type<float>::def;
  V.dimname[2] = axis_name3;
  V.start[2] = axis3(0);
  V.stride[2] = (axis3(V.dimsize[2]-1)-axis3(0))/(V.dimsize[2]-1);
  V.end[2] = axis3(V.dimsize[2]-1);
}

void Diagnostics::savenE(Mapper &mapper, int i,
                         Array1dc &nj, Array1dc &Ej,
                         Array1dc &nk, Array1dc &Ek)
{
  start.resize(2);
  edge.resize(2);

  start[0] = iter.pos(i);
  edge[0] = 1;
  start[1] = 0;
  edge[1] = nj.rows();

  if (njDiag) {
    Array1dr n(nj.size());

    n = real(nj);
    n = mapper.toPhysicalDensity(n);
    hdf.writeSDvar(re_nj_hdf_name, start, edge, n.data());

    n = imag(nj);
    n = mapper.toPhysicalDensity(n);
    hdf.writeSDvar(im_nj_hdf_name, start, edge, n.data());
  }

  if (EjDiag) {
    Array1dr E(Ej.size());

    E = real(Ej);
    E = mapper.toPhysicalElectricField(E);
    hdf.writeSDvar(re_Ej_hdf_name, start, edge, E.data());

    E = imag(Ej);
    E = mapper.toPhysicalElectricField(E);
    hdf.writeSDvar(im_Ej_hdf_name, start, edge, E.data());
  }

  if (nkDiag) {
    int N = mapper.getGridSize();
    // Convert to physical dimension
    Array1dc nkphys( fourier::fftshift(Array1dc( mapper.toPhysicalDensity(nk)*(1.0/N) )) );

    Array1dr n(N);
    n = real(nkphys);
    {
      // Remove zero-padded wavenumbers
      int M = mapper.getZeroPaddingSize();
      Array1dr tmp(1.0*n);
      n.resize(2*M+1);
      n = tmp(blitz::Range(N/2-M,N/2+M));
      edge[1] = n.size();
    }
    hdf.writeSDvar(re_nk_hdf_name, start, edge, n.data());

    n.resize(N);
    n = imag(nkphys);
    {
      // Remove zero-padded wavenumbers
      int M = mapper.getZeroPaddingSize();
      Array1dr tmp(1.0*n);
      n.resize(2*M+1);
      n = tmp(blitz::Range(N/2-M,N/2+M));
      edge[1] = n.size();
    }
    hdf.writeSDvar(im_nk_hdf_name, start, edge, n.data());
  }

  if (EkDiag) {
    int N = mapper.getGridSize();
    Array1dc Ekphys( fourier::fftshift(Array1dc( mapper.toPhysicalElectricField(Ek)*(1.0/N) )) );

    Array1dr E(N);
    E = real(Ekphys);
    {
      // Remove zero-padded wavenumbers
      int M = mapper.getZeroPaddingSize();
      Array1dr tmp(1.0*E);
      E.resize(2*M+1);
      E = tmp(blitz::Range(N/2-M,N/2+M));
      edge[1] = E.size();
    }
    hdf.writeSDvar(re_Ek_hdf_name, start, edge, E.data());

    E.resize(N);
    E = imag(Ekphys);
    {
      // Remove zero-padded wavenumbers
      int M = mapper.getZeroPaddingSize();
      Array1dr tmp(1.0*E);
      E.resize(2*M+1);
      E = tmp(blitz::Range(N/2-M,N/2+M));
      edge[1] = E.size();
    }
    hdf.writeSDvar(im_Ek_hdf_name, start, edge, E.data());
  }
}

void Diagnostics::saveNPH(Mapper &mapper, int i, Real N, Real P, Real H)
{
  if (NPHDiag) {
    start.resize(1);
    edge.resize(1);

    start[0] = iter.pos(i);
    edge[0] = 1;

    hdf.writeSDvar(N_hdf_name, start, edge, &N);
    hdf.writeSDvar(P_hdf_name, start, edge, &P);
    hdf.writeSDvar(H_hdf_name, start, edge, &H);
  }
}

void Diagnostics::saveR(Mapper &mapper, int i, Real R)
{
  if (NPHDiag) {
    start.resize(1);
    edge.resize(1);

    start[0] = iter.pos(i);
    edge[0] = 1;

    hdf.writeSDvar(R_hdf_name, start, edge, &R);
  }
}

void Diagnostics::saveW(Mapper &mapper, int i, Real W)
{
  if (NPHDiag) {
    start.resize(1);
    edge.resize(1);

    start[0] = iter.pos(i);
    edge[0] = 1;

    hdf.writeSDvar(W_hdf_name, start, edge, &W);
  }
}

void Diagnostics::saveDiffusion(Mapper &mapper, int i)
{
  start.resize(2);
  edge.resize(2);

  start[0] = iter.pos(i);
  edge[0] = 1;

  if (FeDiag) {
    Array1dr Fe( mapper.getDiffusionElectronPdf() );
    start[1] = 0;
    edge[1] = Fe.rows();
    hdf.writeSDvar(Fe_hdf_name, start, edge, Fe.data());
  }

  if (nueDiag) {
    Array1dr nue( mapper.getDiffusionLangmuirWavesDamping() );
    start[1] = 0;
    edge[1] = nue.rows();
    hdf.writeSDvar(nue_hdf_name, start, edge, nue.data());
  }
}

void Diagnostics::saveSk2(Mapper &mapper, Spectra &spectra)
{
  if (Sk2Diag && spectra.isKspecDiagnosed()) {
    if (spectra.isKspecIntegrated()) {

      Array1dr nk2, Ek2;
      spectra.getSk2(mapper, nk2, Ek2);

      start.resize(2);
      edge.resize(2);

      start[0] = spectra.getKspecTindex();
      edge[0] = 1;
      start[1] = 0;
      edge[1] = Ek2.numElements();

#if 0
      cout << "start=" << start << " edge=" << edge << endl;
#endif
      hdf.writeSDvar(nk2_hdf_name, start, edge, nk2.data());
      hdf.writeSDvar(Ek2_hdf_name, start, edge, Ek2.data());
    }
  }
}

void Diagnostics::saveSkf2(Mapper &mapper, int i, bool atEnd, Spectra &spectra)
{
  if (Skf2Diag) {
    if (spectra.isIspecDiagnosed()) { // Ion line
      if (spectra.isIspecIntegrated() && atEnd) {

        start.resize(2);
        edge.resize(2);

        Array2dr nkf2;
        spectra.getIkf2(mapper, nkf2);

        start[0] = 0;
        edge[0] = nkf2.rows();
        start[1] = 0;
        edge[1] = nkf2.cols();
#if 0
        cout << "start=" << start << " edge=" << edge << endl;
#endif
        hdf.writeSDvar(nkf2_hdf_name, start, edge, nkf2.data());
      }
      if (! spectra.isIspecIntegrated() && spectra.isIspecComplete(i)) {

        start.resize(3);
        edge.resize(3);

        Array2dr nkf2;
        spectra.getIkf2(mapper, nkf2);

        start[0] = spectra.getIspecTindex();
        edge[0] = 1;
        start[1] = 0;
        edge[1] = nkf2.rows();
        start[2] = 0;
        edge[2] = nkf2.cols();
#if 0
        cout << "start=" << start << " edge=" << edge << endl;
#endif
        hdf.writeSDvar(nkf2_hdf_name, start, edge, nkf2.data());
      }
    }
    if (spectra.isEspecDiagnosed()) { // Electric field spectrum
      if (spectra.isEspecIntegrated() && atEnd) {

        start.resize(2);
        edge.resize(2);

        Array2dr Ekf2;
        spectra.getEkf2(mapper, Ekf2);

        start[0] = 0;
        edge[0] = Ekf2.rows();
        start[1] = 0;
        edge[1] = Ekf2.cols();
#if 0
        cout << "start=" << start << " edge=" << edge << endl;
#endif
        hdf.writeSDvar(Ekf2_hdf_name, start, edge, Ekf2.data());
      }
      if (! spectra.isEspecIntegrated() && spectra.isEspecComplete(i)) {

        start.resize(3);
        edge.resize(3);

        Array2dr Ekf2;
        spectra.getEkf2(mapper, Ekf2);

        start[0] = spectra.getEspecTindex();
        edge[0] = 1;
        start[1] = 0;
        edge[1] = Ekf2.rows();
        start[2] = 0;
        edge[2] = Ekf2.cols();
#if 0
        cout << "start=" << start << " edge=" << edge << endl;
#endif
        hdf.writeSDvar(Ekf2_hdf_name, start, edge, Ekf2.data());
      }
    }
    if (spectra.isUspecDiagnosed()) { // Up plasma line
      if (spectra.isUspecIntegrated() && atEnd) {

        start.resize(2);
        edge.resize(2);

        Array2dr ukf2;
        spectra.getUkf2(mapper, ukf2);

        start[0] = 0;
        edge[0] = ukf2.rows();
        start[1] = 0;
        edge[1] = ukf2.cols();
#if 0
        cout << "start=" << start << " edge=" << edge << endl;
#endif
        hdf.writeSDvar(ukf2_hdf_name, start, edge, ukf2.data());
      }
      if (! spectra.isUspecIntegrated() && spectra.isUspecComplete(i)) {

        start.resize(3);
        edge.resize(3);

        Array2dr ukf2;
        spectra.getUkf2(mapper, ukf2);

        start[0] = spectra.getUspecTindex();
        edge[0] = 1;
        start[1] = 0;
        edge[1] = ukf2.rows();
        start[2] = 0;
        edge[2] = ukf2.cols();
#if 0
        cout << "start=" << start << " edge=" << edge << endl;
#endif
        hdf.writeSDvar(ukf2_hdf_name, start, edge, ukf2.data());
      }
    }
    if (spectra.isDspecDiagnosed()) { // Down plasma line
      if (spectra.isDspecIntegrated() && atEnd) {

        start.resize(2);
        edge.resize(2);

        Array2dr dkf2;
        spectra.getDkf2(mapper, dkf2);

        start[0] = 0;
        edge[0] = dkf2.rows();
        start[1] = 0;
        edge[1] = dkf2.cols();
#if 0
        cout << "start=" << start << " edge=" << edge << endl;
#endif
        hdf.writeSDvar(dkf2_hdf_name, start, edge, dkf2.data());
      }
      if (! spectra.isDspecIntegrated() && spectra.isDspecComplete(i)) {

        start.resize(3);
        edge.resize(3);

        Array2dr dkf2;
        spectra.getDkf2(mapper, dkf2);

        start[0] = spectra.getDspecTindex();
        edge[0] = 1;
        start[1] = 0;
        edge[1] = dkf2.rows();
        start[2] = 0;
        edge[2] = dkf2.cols();
#if 0
        cout << "start=" << start << " edge=" << edge << endl;
#endif
        hdf.writeSDvar(dkf2_hdf_name, start, edge, dkf2.data());
      }
    }
    if (spectra.isSEEspecDiagnosed()) { // SEE spectrum
      if (spectra.isSEEspecIntegrated() && atEnd) {

        start.resize(1);
        edge.resize(1);

        Array2dr seef;
        spectra.getSEEf(mapper, seef);

        start[0] = 0;
        edge[0] = seef.rows();
#if 0
        cout << "start=" << start << " edge=" << edge << endl;
#endif
        hdf.writeSDvar(seef_hdf_name, start, edge, seef.data());
      }
      if (! spectra.isSEEspecIntegrated() && spectra.isSEEspecComplete(i)) {

        start.resize(2);
        edge.resize(2);

        Array2dr seef;
        spectra.getSEEf(mapper, seef);

        start[0] = spectra.getSEEspecTindex();
        edge[0] = 1;
        start[1] = 0;
        // tricky! seef is an 2D array with dim (1,nfft)
        // since first dim contain k=0
        edge[1] = seef.columns();
#if 0
        cout << "start=" << start << " edge=" << edge << endl;
#endif
        hdf.writeSDvar(seef_hdf_name, start, edge, seef.data());
      }
    }
  }
}

#if 0
if (spectra.spec_nt>0)
{
  ofstream id1("nkt.out");
  ofstream id2("Ekt.out");
  for (int it=spectra.spec_start_iter; it<=spectra.spec_end_iter;
       it+=spectra.spec_stride_iter) {
    if (solver.collision) {
      id1 << solver.h*solver.tau*it << '\t';
      id2 << solver.h*solver.tau*it << '\t';
    } else {
      id1 << solver.h*it << '\t';
      id2 << solver.h*it << '\t';
    }
    for (int ik=0; ik<spectra.nb_ks; ++ik) {
      id1 << real(ISnkt(ik, it-spectra.spec_start_iter)) << "\t"
          << imag(ISnkt(ik, it-spectra.spec_start_iter)) << "\t";
      id2 << real(ISEkt(ik, it-spectra.spec_start_iter)) << "\t"
          << imag(ISEkt(ik, it-spectra.spec_start_iter)) << "\t";
    }
    id1 << endl;
    id2 << endl;
  }
}
#endif
