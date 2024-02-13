/**************************************************************************
 *
 * $Id: time-integrate-noFp.cpp,v 1.1 2011/04/08 15:21:38 patrick Exp $
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

void ZakharovSolver::timeIntegrateNoFp()
{
  // from Payne et al., Numerical solution of the Zakharov equations,
  // J. Comp. Phys., 50, 482-498, 1983

  // Eq. 10 Integrate nk from t=0 to t=h without ponderomotive contribution
  nkh(kp) = nk(kp)*sDecay1(kp) + dnk(kp)*sDecay2(kp);
  nkh(kn) = nk(kn)*sDecay1(kn) + dnk(kn)*sDecay2(kn);

  // Add density fluctuation source nk(h) = nk(h) + delta_nk
  calculateNkSource();

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
  nEj = real(nj)*Ej;
  fft.direct(nEj, nEk);
  zeroPadding(nEk);

  // Eq. 13 nEestk integrated from t=0 to t=h
  nEestk(kp) = nEk(kp)*lDecay(kp);
  nEestk(kn) = nEk(kn)*lDecay(kn);
  // Fourier tranform nEestk(h) into nEestj(h)
  fft.inverse(nEestk, nEestj);

  // Compute SSE term
  nEjmean(0) = blitz::sum(nEestj)/Real(N);

  // Eq. 15 Integrate Ej from t=0 to t=h
  Ejh = (Eestj - ihOver2 * nEestj) / (1.0 + ihOver2 * real(njh));
  // Fourier tranform Ej(h) into Ek(h)
  fft.direct(Ejh, Ekh);
  zeroPadding(Ekh);

  // Add electrostatic field fluctuation source Ek(h) = Ek(h) + delta_Ek
  calculateEkSource();

  // Eq. 19 Integrate dnk from t=0 to t=h
  dnk(kp) = dnk(kp)*sDecay4(kp) - hOver2k2(kp)*(nk(kp)*sDecay4(kp) + nkh(kp));
  dnk(kn) = dnk(kn)*sDecay4(kn) - hOver2k2(kn)*(nk(kn)*sDecay4(kn) + nkh(kn));

  // Update nk(h) into nk(0)
  nk(kp) = nkh(kp);
  nk(kn) = nkh(kn);

  // Update Ek(h) into Ek(0)
  Ek(kp) = Ekh(kp);
  Ek(kn) = Ekh(kn);

  // Update nj(h) into nj(0)
  nj = njh;

  // Update dnj(h) into dnj(0)
  zeroPadding(dnk);
  fft.inverse(dnk, dnj);

  // Update Ej(h) into Ej(0)
  zeroPadding(Ek);
  fft.inverse(Ek , Ej);

  if (diffusion) mapper.updateDiffusion(Ek, nue);

  ++iter;
}

