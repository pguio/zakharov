/**************************************************************************
 *
 * $Id: spectra.cpp,v 1.48 2019/05/10 16:54:43 patrick Exp $
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

#include <spectra.h>

using blitz::pow2;

using parser::header;
using parser::yesno;

using std::ostream;
using std::string;

using std::cout;
using std::endl;


ostream& operator<<(ostream& os, const Spectra::kSpectraConf &conf)
{
  if (conf.isDiagnosed())  {
    os
        << "<|S(k)|^2> time iter             = " << conf.iter << '\n'
        << "<|S(k)|^2> wavenumbers index     = " << conf.ik << '\n'
        << "<|S(k)|^2> wavenumbers [m-1]     = " << rangeStatus(conf.k);
  } else {
    os << "no <|S(k)|^2>   diagnostics";
  }
  return os;
}

ostream& operator<<(ostream& os, const Spectra::kfSpectraConf &conf)
{
  const char type[2]= {
    conf.type, '\0'
  };
  const string spec("<|"+string(type)+"(k,f)|^2>");
  if (conf.isDiagnosed()) {

    os
        << spec << " time iter           = " << conf.iter << '\n'

        << spec << " sampling interval   = " << conf.ts << " [s]\n"
        << spec << " sampling freq       = " << conf.fs << " [Hz]\n"
        << spec << " number of pts       = " << conf.ns << '\n'

        << spec << " wave numbers index  = " << conf.ik << '\n'
        << spec << " wave numbers (m-1)  = " << rangeStatus(conf.k) << '\n'

        << spec << " integrate           = " << yesno(conf.integrate) << '\n'
        << spec << " post-integration nb = " << conf.postint <<'\n'

        << spec << " FFT setup           = " << conf.fft;
  } else {
    os << "no " << spec << " diagnostics";
  }
  return os;
}

ostream& operator<<(ostream& os, const Spectra &spectra)
{
  return os
         << header("Spectra setup")
         << spectra.kspec << "\n##\n"
         << spectra.ispec << "\n##\n"
         << spectra.Espec << "\n##\n"
         << spectra.uspec << "\n##\n"
         << spectra.dspec << "\n##\n"
         << spectra.seespec;
}


Spectra::Spectra(int nargs, char *args[]) : Parser(nargs, args),
  ispec('n'), uspec('u'), dspec('d'), Espec('E'), seespec('s')
{
  initParsing(nargs, args);
  parseParam();
}

Spectra::~Spectra()
{}


void Spectra::initParsing(int nargs, char *args[])
{
  const string id(PACKAGE);
  const string version(
    VERSION "\n"
    "$Id: spectra.cpp,v 1.48 2019/05/10 16:54:43 patrick Exp $"
    "\n");
  const string copyright(ZAK_COPYRIGHT);

  registerClass("Spectra");
  registerPackage(id, version, copyright);

  using parser::types::boolean;
  using parser::types::integer;
  using parser::types::intVect;

  // wavenumber spectra
  insertOption(_kspec_iter, "kspec_iter", intVect, "<|S(k)|^2> time start,stride,end", Any(kspec.iter));
  insertOption(_kspec_ik  , "kspec_ik"  , intVect, "<|S(k)|^2> k start,stride,end"   , Any(kspec.ik));

  // ion line
  insertOption(_ispec_iter, "ispec_iter", intVect, "<|n(k,f)|^2> time start,stride,end", Any(ispec.iter));
  insertOption(_ispec_ik  , "ispec_ik"  , intVect, "<|n(k,f)|^2> k start,stride,end"   , Any(ispec.ik));
  insertOption(_ispec_ns  , "ispec_ns"  , integer, "Number of freq points <|n(k,f)|^2>", Any(ispec.ns));
  insertOption(_ispec_integrate, "ispec_integrate", boolean, "Integrate <|n(k,f)|^2>", Any(ispec.integrate));
  insertOption(_ispec_postint  , "ispec_postint"  , integer, "Post Integration number <|n(k,f)|^2>", Any(ispec.postint));

  // electric field spectra
  insertOption(_Espec_iter, "Espec_iter", intVect, "<|E(k,f)|^2> time start,stride,end", Any(Espec.iter));
  insertOption(_Espec_ik  , "Espec_ik"  , intVect, "<|E(k,f)|^2> k start,stride,end"   , Any(Espec.ik));
  insertOption(_Espec_ns  , "Espec_ns"  , integer, "Number of freq points <|E(k,f)|^2>", Any(Espec.ns));
  insertOption(_Espec_integrate, "Espec_integrate", boolean, "Integrate <|E(k,f)|^2>", Any(Espec.integrate));
  insertOption(_Espec_postint  , "Espec_postint"  , integer, "Post Integration number <|E(k,f)|^2>", Any(Espec.postint));

  // up plasma line
  insertOption(_uspec_iter, "uspec_iter", intVect, "<|u(k,f)|^2> time start,stride,end", Any(uspec.iter));
  insertOption(_uspec_ik  , "uspec_ik"  , intVect, "<|u(k,f)|^2> k start,stride,end", Any(uspec.ik));
  insertOption(_uspec_ns  , "uspec_ns"  , integer, "Number of freq points <|u(k,f)|^2>",Any(uspec.ns));
  insertOption(_uspec_integrate, "uspec_integrate", boolean, "Integrate <|u(k,f)|^2>", Any(uspec.integrate));
  insertOption(_uspec_postint, "uspec_postint", integer, "Post Integration number <|u(k,f)|^2>", Any(uspec.postint));

  // down plasma line
  insertOption(_dspec_iter, "dspec_iter", intVect, "<|d(k,f)|^2> time start,stride,end", Any(dspec.iter));
  insertOption(_dspec_ik, "dspec_ik", intVect, "<|d(k,f)|^2> k start,stride,end", Any(dspec.ik));
  insertOption(_dspec_ns, "dspec_ns", intVect, "Number of freq points <|d(k,f)|^2>", Any(dspec.ns));
  insertOption(_dspec_integrate, "dspec_integrate", boolean, "Integrate <|d(k,f)|^2>", Any(dspec.integrate));
  insertOption(_dspec_postint, "dspec_postint", integer, "Post Integration number <|d(k,f)|^2>", Any(dspec.postint));

  // SEE spectrum
  insertOption(_seespec_iter, "seespec_iter", intVect, "<|SEE(f)|> time start,stride,end", Any(seespec.iter));
  insertOption(_seespec_ns, "seespec_ns", intVect, "Number of freq points <|SEE(f)|>", Any(seespec.ns));
  insertOption(_seespec_integrate, "seespec_integrate", boolean, "Integrate <|SEE(f)|>", Any(seespec.integrate));
  insertOption(_seespec_postint, "seespec_postint", integer, "Post Integration number <|SEE(f)|>", Any(seespec.postint));

}

void Spectra::parseParam()
{
  parseOption(_kspec_iter, kspec.iter.data);
  parseOption(_kspec_ik, kspec.ik.data);

  parseOption(_ispec_iter, ispec.iter.data);
  parseOption(_ispec_ik, ispec.ik.data);
  parseOption(_ispec_ns, ispec.ns);
  parseOption(_ispec_integrate, ispec.integrate);
  parseOption(_ispec_postint, ispec.postint);

  parseOption(_Espec_iter, Espec.iter.data);
  parseOption(_Espec_ik, Espec.ik.data);
  parseOption(_Espec_ns, Espec.ns);
  parseOption(_Espec_integrate, Espec.integrate);
  parseOption(_Espec_postint, Espec.postint);

  parseOption(_uspec_iter, uspec.iter.data);
  parseOption(_uspec_ik, uspec.ik.data);
  parseOption(_uspec_ns, uspec.ns);
  parseOption(_uspec_integrate, uspec.integrate);
  parseOption(_uspec_postint, uspec.postint);

  parseOption(_dspec_iter, dspec.iter.data);
  parseOption(_dspec_ik, dspec.ik.data);
  parseOption(_dspec_ns, dspec.ns);
  parseOption(_dspec_integrate, dspec.integrate);
  parseOption(_dspec_postint, dspec.postint);

  parseOption(_seespec_iter, seespec.iter.data);
  parseOption(_seespec_ns, seespec.ns);
  parseOption(_seespec_integrate, seespec.integrate);
  parseOption(_seespec_postint, seespec.postint);
}

void Spectra::initialise(Mapper &mapper)
{
  if (kspec.isDiagnosed())
    kspec.initialise(mapper);

  if (ispec.isDiagnosed())
    ispec.initialise(mapper);

  if (Espec.isDiagnosed())
    Espec.initialise(mapper);

  if (uspec.isDiagnosed())
    uspec.initialise(mapper);

  if (dspec.isDiagnosed())
    dspec.initialise(mapper);

  if (seespec.isDiagnosed()) {
    // find zero-centered index of k=0
    int ik0 = mapper.getGridSize()/2;
    seespec.ik = RangeSpec<int>(ik0, 1, ik0);
    seespec.initialise(mapper);
  }
}

void Spectra::updateSpectra(Mapper &mapper, int i,
                            const Array1dc &nk, const Array1dc &Ek,
                            const Array1dc &nE0)
{
  if (kspec.inSpectra(i))
    kspec.update(i, nk, Ek);

  if (ispec.inSpectra(i))
    ispec.update(i, nk);

  if (Espec.inSpectra(i))
    Espec.update(i, Ek);

  if (uspec.inSpectra(i))
    uspec.update(i, Ek);

  if (dspec.inSpectra(i))
    dspec.update(i, Ek);

  if (nE0.size()>0 && seespec.inSpectra(i))
    seespec.update(i, nE0);
}

void Spectra::getSk2(Mapper &mapper, Array1dr &nk2, Array1dr &Ek2) const
{
  // Parseval'theorem  \sum |E_j|^2 = 1/nfft \sum |E_k|^2
  int nk = mapper.getGridSize();

  nk2.resize(kspec.nk2.size());
  nk2 = kspec.nk2*(1.0/pow2(nk));

  if (kspec.nspec > 1) nk2 = nk2*(1.0/kspec.nspec);

  Ek2.resize(kspec.Ek2.size());
  Ek2 = kspec.Ek2*(1.0/pow2(nk));

  if (kspec.nspec > 1) Ek2 = Ek2*(1.0/kspec.nspec);

  // Normalise to physical quantities
  nk2 = mapper.toPhysicalDensity2(nk2);
  Ek2 = mapper.toPhysicalElectricField2(Ek2);
}

void Spectra::getIkf2(Mapper &mapper, Array2dr &nkf2) const
{
  int nk = mapper.getGridSize();

  nkf2.resize(ispec.skf2.shape());
  nkf2 = ispec.skf2*(1.0/pow2(nk)/pow2(ispec.ns));
  nkf2 = mapper.toPhysicalDensity2(nkf2);

  if (ispec.nspec > 1) nkf2 = nkf2*(1.0/ispec.nspec);

  // center zero-frequency
  for (unsigned ik=0; ik<ispec.k.size(); ++ik) {
    Array1dr skf2(ispec.ns);
    skf2 = nkf2(ik,Range::all());
    skf2 = fourier::fftshift(skf2);
    nkf2(ik,Range::all()) = skf2;
  }
}

void Spectra::getEkf2(Mapper &mapper, Array2dr &Ekf2) const
{
  int nk = mapper.getGridSize();

  Ekf2.resize(Espec.skf2.shape());
  Ekf2 = Espec.skf2*(1.0/pow2(nk)/pow2(Espec.ns));
  Ekf2 = mapper.toPhysicalElectricField2(Ekf2);

  if (Espec.nspec > 1) Ekf2 *= (1.0/Espec.nspec);

  // center zero-frequency
  for (unsigned ik=0; ik<Espec.k.size(); ++ik) {
    Array1dr skf2(Espec.ns);
    skf2 = Ekf2(ik,Range::all());
    skf2 = fourier::fftshift(skf2);
    Ekf2(ik,Range::all()) = skf2;
  }
}

void Spectra::getUkf2(Mapper &mapper, Array2dr &ukf2) const
{
  int nk = mapper.getGridSize();

  ukf2.resize(uspec.skf2.shape());
  ukf2 = uspec.skf2*(1.0/pow2(nk)/pow2(uspec.ns));
  ukf2 = mapper.toPhysicalElectricField2(ukf2);

  if (uspec.nspec > 1) ukf2 = ukf2*(1.0/uspec.nspec);

  // center zero-frequency and transform Ek2 into k^2 Ek2
  for (unsigned ik=0; ik<uspec.k.size(); ++ik) {
    Array1dr skf2(uspec.ns);
    skf2 = ukf2(ik,Range::all());
    skf2 = fourier::fftshift(skf2);
    skf2 = pow2(uspec.k(ik))*skf2;
    ukf2(ik,Range::all()) = skf2;
  }
  ukf2 = mapper.toPhysicalPlasmaLine(ukf2);
}

void Spectra::getDkf2(Mapper &mapper, Array2dr &dkf2) const
{
  int nk = mapper.getGridSize();

  dkf2.resize(dspec.skf2.shape());
  dkf2 = dspec.skf2*(1.0/pow2(nk)/pow2(dspec.ns));
  dkf2 = mapper.toPhysicalElectricField2(dkf2);

  if (dspec.nspec > 1) dkf2 = dkf2*(1.0/dspec.nspec);

  // center zero-frequency and transform Ek2 into k^2 Ek2
  for (unsigned ik=0; ik<dspec.k.size(); ++ik) {
    Array1dr skf2(dspec.ns);
    skf2 = dkf2(ik,Range::all());
    skf2 = fourier::fftshift(skf2);
    skf2 = pow2(dspec.k(ik))*skf2;
    dkf2(ik,Range::all()) = skf2;
  }
  dkf2 = mapper.toPhysicalPlasmaLine(dkf2);
}

void Spectra::getSEEf(Mapper &mapper, Array2dr &nEf) const
{
  // Parseval'theorem  \sum n_jE_j^* = 1/nfft \sum n_kE_k^*
  int nk = mapper.getGridSize();

  nEf.resize(seespec.skf2.shape());
  nEf = seespec.skf2*(1.0/pow2(nk)/pow2(seespec.ns));
  nEf = mapper.toPhysicalElectricField2(nEf);

  if (seespec.nspec > 1) nEf = nEf*(1.0/seespec.nspec);

  // center zero-frequency and transform Ek2 into k^2 Ek2
  for (unsigned ik=0; ik<seespec.k.size(); ++ik) {
    Array1dr skf2(seespec.ns);
    skf2 = nEf(ik,Range::all());
    skf2 = fourier::fftshift(skf2);
    skf2 = sqrt(skf2);
    nEf(ik,Range::all()) = skf2;
  }
}
