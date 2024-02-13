/**************************************************************************
 *
 * $Id: time-integrate-fixed-E.cpp,v 1.2 2011/03/26 07:47:28 patrick Exp $
 *
 * Copyright (c) 2011 Patrick Guio <patrick.guio@gmail.com>
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

void ZakharovSolver::timeIntegrateFixedE()
{
  // from Payne et al., Numerical solution of the Zakharov equations,
  // J. Comp. Phys., 50, 482-498, 1983

  // Term (|E(0)|^2)_k in Eq. 10
  Ej2 = zip(pow2(real(Ej))+pow2(imag(Ej)), 0.0, Complex());
  fft.direct(Ej2, Ek2);
  zeroPadding(Ek2);

  // Eq. 10 Integrate nk from t=0 to t=h
  nkh(kp) = nk(kp)*sDecay1(kp) + dnk(kp)*sDecay2(kp) - Ek2(kp)*sDecay3(kp);
  nkh(kn) = nk(kn)*sDecay1(kn) + dnk(kn)*sDecay2(kn) - Ek2(kn)*sDecay3(kn);

#if 1
  // Add density fluctuation source nk(h) = nk(h) + delta_nk
  calculateNkSource();
#endif

  // Fourier tranform nk(h) into nj(h) needed for Eq. 15
  zeroPadding(nkh);
  fft.inverse(nkh, njh);

#if 1
  // Eq. 19 Integrate dnk from t=0 to t=h
  dnk(kp) = dnk(kp)*sDecay4(kp) - hOver2k2(kp)*((nk(kp)+Ek2(kp))*sDecay4(kp) + nkh(kp)+Ek2(kp));
  dnk(kn) = dnk(kn)*sDecay4(kn) - hOver2k2(kn)*((nk(kn)+Ek2(kn))*sDecay4(kn) + nkh(kn)+Ek2(kn));
#else
  dnk(kp) = 0.0;
  dnk(kn) = 0.0;
#endif

  // Update nk(h) into nk(0)
  nk(kp) = nkh(kp);
  nk(kn) = nkh(kn);

  // Update nj(h) into nj(0)
  zeroPadding(nk);
  fft.inverse(nk , nj);
#if 1
  // Update dnj(h) into dnj(0)
  zeroPadding(dnk);
  fft.inverse(dnk, dnj);
#endif

  ++iter;
}

