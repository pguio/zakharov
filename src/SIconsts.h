/**************************************************************************
 *
 * $Id: SIconsts.h,v 1.6 2011/03/26 07:47:28 patrick Exp $
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

#ifndef SICONSTS_H
#define SICONSTS_H

class SI {
public:

  // Universal constants
  static const double c;      // m s-1 speed of light
  static const double G;      // m^3 kg-1 s-2 gravitation constant
  static const double h;      // J s      Planck constant
  static const double hbar;   // J s      Planck constant/2pi

  // Electromagnetic constants
  static const double e;     // C elementary charge
  static const double mu0;   // N A-2    vacuum permeability
  static const double eps0;  // F m-1 vacuum permitivity

  // Atomic constants
  static const double me;    // kg electron mass
  static const double mp;    // kg      proton mass
  static const double mn;    // kg      neutron mass

  // Thermodynamic constants
  static const double Kb;    // J K-1 Boltzmann constant
  static const double L;     // mol-1    Avogadro constant
  static const double R;     // J mol-1 K-1    Gas constant

  // Conversion factors
  static const double eV;    // J electron volt
  static const double amu;   // kg atomic mass unit
  static const double atm;   // Pa    standard athmoshere
  static const double cal;   // J     thermodynamic calorie
};

#endif // SICONSTS_H
