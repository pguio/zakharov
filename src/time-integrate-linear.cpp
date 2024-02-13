/**************************************************************************
 *
 * $Id: time-integrate-linear.cpp,v 1.3 2011/04/08 15:29:02 patrick Exp $
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


void ZakharovSolver::timeIntegrateLinear()
{
  // from Payne et al., Numerical solution of the Zakharov equations,
  // J. Comp. Phys., 50, 482-498, 1983

  // Eq. 10 Integrate nk from t=0 to t=h
  // The last (nonlinear) term for ponderomotive force is omitted
  nkh(kp) = nk(kp)*sDecay1(kp) + dnk(kp)*sDecay2(kp);
  nkh(kn) = nk(kn)*sDecay1(kn) + dnk(kn)*sDecay2(kn);

  // Add density fluctuation source nk(h) = nk(h) + delta_nk
  calculateNkSource();

  if (diffusion) calculateLangmuirDecay();

  // Eq. 9 Integrate Ek from t=0 to t=h
  // The last (nonlinear) term for density inhomogeneity is omitted
  Ek(kp) *= lDecay(kp);
  Ek(kn) *= lDecay(kn);

  // Add electrostatic field fluctuation source Ek(h) = Ek(h) + delta_Ek
  calculateEkSource();

  // Eq. 19 Integrate dnk from t=0 to t=h
  // The (nonlinear) terms for ponderomotive force in the last two terms
  // are omitted
  dnk(kp) = dnk(kp)*sDecay4(kp)-hOver2k2(kp)*(nk(kp)*sDecay4(kp)+nkh(kp));
  dnk(kn) = dnk(kn)*sDecay4(kn)-hOver2k2(kn)*(nk(kn)*sDecay4(kn)+nkh(kn));

  // Update nk(h) into nk(0)
  nk(kp) = nkh(kp);
  nk(kn) = nkh(kn);

  // Update nj(h) into nj(0)
  zeroPadding(nk);
  fft.inverse(nk , nj);

  // Update dnj(h) into dnj(0)
  zeroPadding(dnk);
  fft.inverse(dnk, dnj);

  // Update Ej(h) into Ej(0)
  zeroPadding(Ek);
  fft.inverse(Ek , Ej);

  // Compute SSE term
  nEj = real(nj)*Ej;
  nEjmean(0) = blitz::sum(nEj)/Real(N);

  if (diffusion) mapper.updateDiffusion(Ek, nue);

  ++iter;
}

