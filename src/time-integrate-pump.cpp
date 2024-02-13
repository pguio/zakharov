/**************************************************************************
 *
 * $Id: time-integrate-pump.cpp,v 1.3 2011/11/08 15:14:24 patrick Exp $
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

#include <zakharov.h>

using blitz::imag;
using blitz::pow2;
using blitz::real;
using blitz::zip;

using std::cout;
using std::endl;

void ZakharovSolver::timeIntegratePump()
{

  // Adapted from Payne et al., Numerical solution of the Zakharov equations,
  // J. Comp. Phys., 50, 482-498, 1983
  // and modified to take into account a pump drive Ep at k=0

  // Pump values at time t and t+h
  Ep = mapper.getPump(h*iter);
  Eph = mapper.getPump(h*(iter+1));

  // Term (|E(0)|^2)_k in Eq. 10
  Ej2 = zip(pow2(blitz::real(Ej+Ep))+pow2(blitz::imag(Ej+Ep)),0.0,Complex());
  fft.direct(Ej2, Ek2);
  zeroPadding(Ek2);

  // Eq. 10 Integrate nk from t=0 to t=h
  nkh(kp) = nk(kp)*sDecay1(kp) + dnk(kp)*sDecay2(kp) - Ek2(kp)*sDecay3(kp);
  nkh(kn) = nk(kn)*sDecay1(kn) + dnk(kn)*sDecay2(kn) - Ek2(kn)*sDecay3(kn);

#if 1
  // Add density fluctuation source  nk(h) = nk(h) + delta_nk
  calculateNkSource();
#endif

  // Fourier tranform nk(h) into nj(h) needed for Eq. 15
  zeroPadding(nkh);
  fft.inverse(nkh, njh);


  if (diffusion) calculateLangmuirDecay();

  // Eq. 12 Eestk integrated from t=0 to t=h
  Eestk(kp) = Ek(kp)*lDecay(kp);
  Eestk(kn) = Ek(kn)*lDecay(kn);
  // Fourier tranform Eestk(h) into Eestj(h)
  zeroPadding(Eestk);
  fft.inverse(Eestk, Eestj);

  // Term (n(0)E(0))_k in Eq. 13
  nEj = real(nj)*(Ej + Ep);
  fft.direct(nEj, nEk);
  zeroPadding(nEk);

  // Eq. 13
  nEestk(kp) = nEk(kp)*lDecay(kp);
  nEestk(kn) = nEk(kn)*lDecay(kn);
  // Fourier tranform nEk into nEj
  zeroPadding(nEestk);
  fft.inverse(nEestk, nEestj);
#if 1
  nEjmean(0) = blitz::sum(nEestj)/Real(N);
  nEestj -= nEjmean(0);
#endif

  // Eq. 15 Integrate Ej from t=0 to t=h
  Ejh = (Eestj - ihOver2*(nEestj+real(njh)*Eph)) / (1.0 + ihOver2*real(njh));
  // Fourier tranform Ej(h) into Ek(h)
  fft.direct(Ejh, Ekh);
  zeroPadding(Ekh);
  Ekh(0) = 0.0;

#if 1
  // Add electrostatic field fluctuation source Ek(h) = Ek(h) + delta_Ek
  calculateEkSource();
#endif

  // Term (|E(h)||^2)_k in Eq. 19
  Ejh2 = zip(pow2(real(Ejh))+pow2(imag(Ejh)), 0.0, Complex());
  fft.direct(Ejh2, Ekh2);
  zeroPadding(Ekh2);
#if 0
  Ekh2(0) = 0.0;
#endif

#if 1
  // Eq. 19 Integrate dnk from t=0 to t=h
  dnk(kp) = dnk(kp)*sDecay4(kp) - hOver2k2(kp)*((nk(kp)+Ek2(kp))*sDecay4(kp) + nkh(kp)+Ekh2(kp));
  dnk(kn) = dnk(kn)*sDecay4(kn) - hOver2k2(kn)*((nk(kn)+Ek2(kn))*sDecay4(kn) + nkh(kn)+Ekh2(kn));
#else
  dnk(kp) = 0.0;
  dnk(kn) = 0.0;
#endif

  nk(kp) = nkh(kp);
  nk(kn) = nkh(kn);

  Ek(kp) = Ekh(kp);
  Ek(kn) = Ekh(kn);

  zeroPadding(nk);
  fft.inverse(nk , nj);
#if 1
  zeroPadding(dnk);
  fft.inverse(dnk, dnj);
#endif
  zeroPadding(Ek);
  fft.inverse(Ek , Ej);

  if (diffusion) mapper.updateDiffusion(Ek, nue);

  ++iter;
}
