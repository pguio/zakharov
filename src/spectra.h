/**************************************************************************
 *
 * $Id: spectra.h,v 1.68 2019/05/10 16:54:43 patrick Exp $
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

#ifndef SPECTRA_H
#define SPECTRA_H

#include <fstream>
#include <iostream>

#include <zakharov-defs.h>
#include <mapper.h>

class Spectra : public parser::Parser {
public:

  typedef std::ostream ostream;

  typedef blitz::Range Range;

  typedef fourier::ODFT1D<fourier::complex,fourier::complex> FT1d;


private:

  enum parser_enum { _kspec_iter=1, _kspec_ik, // wavenumber spectra
                     _ispec_iter, _ispec_ik, _ispec_ns, // IL spectra
                     _ispec_integrate, _ispec_postint,

                     _Espec_iter, _Espec_ik, _Espec_ns, // E-field spectra
                     _Espec_integrate, _Espec_postint,

                     _uspec_iter, _uspec_ik, _uspec_ns, // UPL spectra
                     _uspec_integrate, _uspec_postint,

                     _dspec_iter, _dspec_ik, _dspec_ns, // DPL spectra
                     _dspec_integrate, _dspec_postint,

                     _seespec_iter, _seespec_ik, _seespec_ns, // SEE spectra
                     _seespec_integrate, _seespec_postint
                   };

  struct kSpectraConf { // <|S(k)|^2> variables
    RangeSpec<int> iter;
    RangeSpec<int> ik;
    Range ikRange;
    Array1dr k;
    Array1dr t;
    Array1dr nk2;
    Array1dr Ek2;
    int nspec;
    int ispec;
    kSpectraConf() : iter(-1), ik(-1), nspec(0), ispec(-1) {}

    void initialise(Mapper &mapper) {
      using blitz::tensor::i;

      ikRange.setRange(ik.start(),ik.end(),ik.stride());
      k.resize(ikRange.length());
      Array1dr allk ( mapper.getPhysicalWaveNumbers() );
      k = allk(ikRange);

      t.resize(iter.size()-1); // last element of the range not included
      t = (iter.start()+i*iter.stride())*mapper.getPhysicalTimeStep();

      nk2.resize(k.size());
      nk2 = 0.0;
      Ek2.resize(k.size());
      Ek2 = 0.0;
    }
    bool isDiagnosed() const {
      return iter != -1;
    }
    bool inSpectra(int i) const {
      return iter.isWithin(i);
    }
    void update(int i, const Array1dc &nk, const Array1dc &Ek) {
      using blitz::imag;
      using blitz::pow2;
      using blitz::real;

      using std::cout;
      using std::endl;

#if 0
      cout <<  "  update(" << i << ")";
#endif
      int is = (i-iter.start()) % iter.stride();
#if 0
      cout << " is=" << is;
#endif

      if (is == 0) { // Reset arrays and counter to zero
        nk2 = 0.0;
        Ek2 = 0.0;
        nspec = 0;
      }

      Array1dr tmp(nk.size());
      // Accumulate |n(k)|^2 with centered k=0
      tmp = fourier::fftshift(Array1dr( pow2(real(nk))+pow2(imag(nk)) ) );
      nk2 += tmp(ikRange);

      // Accumulate |E(k)|^2 with centered k=0
      tmp = fourier::fftshift(Array1dr( pow2(real(Ek))+pow2(imag(Ek)) ) );
      Ek2 += tmp(ikRange);

      ++nspec;

      // if integrated
      if (is == iter.stride()-1) {
        ++ispec;
      }
#if 0
      cout << " ispec=" << ispec << " nspec=" << nspec << endl;
#endif
    }
  };
  friend ostream& operator<<(ostream& os,
                             const Spectra::kSpectraConf &conf);

  struct kfSpectraConf { // <|S(k,f)|^2> variables
    const char type;
    RangeSpec<int> iter;
    RangeSpec<int> ik;
    Range ikRange;
    int ns;
    bool integrate;
    int postint;
    Array1dr k;
    Array1dr f;
    Array1dr t;
    Array3dc skt;
    Array2dr skf2;
    int nspec;
    int ispec;
    Real ts;
    Real fs;
    FT1d fft;
    kfSpectraConf(const char _type) : type(_type),
      iter(-1), ik(-1), ns(1), integrate(true), postint(1),
      nspec(0), ispec(-1) {}

    void initialise(Mapper &mapper) {
      using blitz::tensor::i;

      ikRange.setRange(ik.start(),ik.end(),ik.stride());
      k.resize(ikRange.length());
      Array1dr allk ( mapper.getPhysicalWaveNumbers() );
      k = allk(ikRange);

      ts = iter.stride()*mapper.getPhysicalTimeStep();
      fs = 1.0/ts;
      f.resize(ns);
      f = fs*(i-ns/2)/ns;

      if (! integrate) {
        int stride = ns*iter.stride()*postint;
        if ((iter.end()-iter.start()) < stride) {
          ostringstream os, classname;
          os << "Problem with time sampling specs:\n"
             << "  Calculated stride too large compared to timespan:\n"
             << "  ns*stride*postint = " << stride
             << ", end-start = " << iter.end()-iter.start();
          classname << "kfSpectraConf<" << type << ">";
          throw ClassException(classname.str(), os.str());
        }
        t.resize((iter.end()-iter.start())/stride);
        t = (iter.start()+i*stride)*mapper.getPhysicalTimeStep();
      }

      skt.resize(k.size(), f.size(), iter.stride());
      skf2.resize(k.size(), f.size());
      skf2 = 0.0;
      if (type == 'n' || type == 'E' || type == 'u' || type == 's')
        fft.resize(ns, 1);
      if (type == 'd')
        fft.resize(ns, -1);
    }
    bool isDiagnosed() const {
      return iter != -1;
    }
    bool inSpectra(int i) const {
      return iter.isWithin(i);
    }
    bool isSpectraComplete() const {
      return nspec == iter.stride()*postint;
    }
    void update(int i, const Array1dc &sk) {
      using blitz::imag;
      using blitz::pow2;
      using blitz::real;

      using std::cout;
      using std::endl;

      using fourier::binom4filter;
#if 0
      cout << type << " update(" << i << ")";
#endif
      if (ns == -1) {
        // TODO: Save time series!
        Array1dc tmp( fourier::fftshift(sk) );

        // Accumulate time series of s(k) with centered k=0
        Array1dc skTimeSeries(ikRange.length());
        skTimeSeries = tmp(ikRange);

      } else {
        // sub-spectra index
        int is = (i-iter.start()) % iter.stride();
#if 0
        cout << " is=" << is;
#endif
        int it = ((i-iter.start()) / iter.stride()) % ns;
#if 0
        cout << " it=" << it;
#endif
        if (type == 's') {
          // dirty trick for SEE since sk.size() = 1
          skt(Range::all(), it, is) = sk(0);
        } else {
          // Accumulate s(k,t) with centered k=0
          Array1dc tmp( fourier::fftshift(sk) );
          skt(Range::all(), it, is) = tmp(ikRange);
        }
        if (!integrate && this->isSpectraComplete() && is == 0 && it == 0) {
          skf2 = 0.0;
          nspec = 0;
        }
        // if time sequence completed
        if (it == ns-1) {
          Array1dc st(ns);
          Array1dc sf(ns);
          for (unsigned ikk=0; ikk<k.size(); ++ikk) {
            // s(k,t) to |s(k,f)|^2 and accumulate it
            st = skt(ikk, Range::all(), is);
            st = binom4filter(st, Complex(0.0,0.0));
            fft.direct(st, sf);
            skf2(ikk, Range::all()) += pow2(real(sf))+pow2(imag(sf));
          }
          ++nspec;
          if (!integrate && this->isSpectraComplete() && is == iter.stride()-1)
            ++ispec;
        }
#if 0
        cout << " ispec=" << ispec << " nspec=" << nspec << endl;
#endif
      }
    }
  };
  friend ostream& operator<<(ostream& os,
                             const Spectra::kfSpectraConf &conf);

public:
  Spectra(int nargs, char *args[]);
  ~Spectra();

  void initialise(Mapper &mapper);

  bool isKspecDiagnosed() const {
    return kspec.isDiagnosed();
  }
  bool isIspecDiagnosed() const {
    return ispec.isDiagnosed();
  }
  bool isEspecDiagnosed() const {
    return Espec.isDiagnosed();
  }
  bool isUspecDiagnosed() const {
    return uspec.isDiagnosed();
  }
  bool isDspecDiagnosed() const {
    return dspec.isDiagnosed();
  }
  bool isSEEspecDiagnosed() const {
    return seespec.isDiagnosed();
  }

  bool isKspecIntegrated() const {
    return kspec.nspec == kspec.iter.stride();
  }
  bool isIspecIntegrated() const {
    return ispec.integrate;
  }
  bool isEspecIntegrated() const {
    return Espec.integrate;
  }
  bool isUspecIntegrated() const {
    return uspec.integrate;
  }
  bool isDspecIntegrated() const {
    return dspec.integrate;
  }
  bool isSEEspecIntegrated() const {
    return seespec.integrate;
  }


  bool isIspecComplete(int i) const {
    return ispec.isSpectraComplete();
  }
  bool isEspecComplete(int i) const {
    return Espec.isSpectraComplete();
  }
  bool isUspecComplete(int i) const {
    return uspec.isSpectraComplete();
  }
  bool isDspecComplete(int i) const {
    return dspec.isSpectraComplete();
  }
  bool isSEEspecComplete(int i) const {
    return seespec.isSpectraComplete();
  }


  int getKspecTindex() const {
    return kspec.ispec;
  }
  int getIspecTindex() const {
    return ispec.ispec;
  }
  int getEspecTindex() const {
    return Espec.ispec;
  }
  int getUspecTindex() const {
    return uspec.ispec;
  }
  int getDspecTindex() const {
    return dspec.ispec;
  }
  int getSEEspecTindex() const {
    return seespec.ispec;
  }


  Array1dr getKspecK() const {
    return kspec.k;
  }
  Array1dr getIspecK() const {
    return ispec.k;
  }
  Array1dr getEspecK() const {
    return Espec.k;
  }
  Array1dr getUspecK() const {
    return uspec.k;
  }
  Array1dr getDspecK() const {
    return dspec.k;
  }
  Array1dr getSEEspecK() const {
    return seespec.k;
  }


  Array1dr getIspecF() const {
    return ispec.f;
  }
  Array1dr getEspecF() const {
    return Espec.f;
  }
  Array1dr getUspecF() const {
    return uspec.f;
  }
  Array1dr getDspecF() const {
    return dspec.f;
  }
  Array1dr getSEEspecF() const {
    return seespec.f;
  }


  Array1dr getKspecT() const {
    return kspec.t;
  }
  Array1dr getIspecT() const {
    return ispec.t;
  }
  Array1dr getEspecT() const {
    return Espec.t;
  }
  Array1dr getUspecT() const {
    return uspec.t;
  }
  Array1dr getDspecT() const {
    return dspec.t;
  }
  Array1dr getSEEspecT() const {
    return seespec.t;
  }


  void updateSpectra(Mapper &mapper, int i,
                     const Array1dc &nk, const Array1dc &Ek,
                     const Array1dc &nE0 = Array1dc());

  void getSk2(Mapper &mapper, Array1dr &nk2, Array1dr &Ek2) const;
  void getIkf2(Mapper &mapper, Array2dr &nkf2) const;
  void getEkf2(Mapper &mapper, Array2dr &Ekf2) const;
  void getUkf2(Mapper &mapper, Array2dr &ukf2) const;
  void getDkf2(Mapper &mapper, Array2dr &dkf2) const;
  void getSEEf(Mapper &mapper, Array2dr &nEf) const;

  friend ostream& operator<<(ostream& os, const Spectra &spectra);

protected:

  typedef parser::Parser Parser;

private:
  kSpectraConf kspec;                // n and E k-spectrum
  kfSpectraConf ispec, uspec, dspec; // ion line, up and down plasma line
  kfSpectraConf Espec;               // electric field spectra
  kfSpectraConf seespec;             // SEE spectra

  void initParsing(int nargs, char *args[]);
  void parseParam();
};

struct rangeStatus {
  typedef std::ostream ostream;
  int n;
  Real start;
  Real end;
  Real stride;
  rangeStatus(const Array1dr &range) :
    n(range.size()), start(range(0)), end(range(range.size()-1)), stride(0)  {
    if (n>1)
      stride = (range(n-1)-range(0))/Real(n-1);
  }
  friend ostream& operator<<(ostream &os, const rangeStatus &val) {
    int p = os.precision();
    os.precision(4);
    if (val.n > 1 ) {
      os << val.start << ":" << val.stride << ":" << val.end
         << std::setprecision(p);
    } else {
      os << val.start << std::setprecision(p);
    }
    return os;
  }
};

#endif // SPECTRA_H

