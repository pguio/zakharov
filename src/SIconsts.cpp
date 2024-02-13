/**************************************************************************
 *
 * $Id: SIconsts.cpp,v 1.10 2011/03/26 07:47:28 patrick Exp $
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

#include <SIconsts.h>

const double m_pi = 3.14159265358979323846;

// Universal constants
const double SI::c    = 299792458.0;     // m s-1        speed of light
const double SI::G    = 6.67259e-11;     // m^3 kg-1 s-2 gravitation constant
const double SI::h    = 6.6260755e-34;   // J s          Planck constant
const double SI::hbar = SI::h/m_pi/2.0;  // J s          Planck constant/2pi

// Electromagnetic constants
const double SI::e    = 1.60217733e-19;             // C     elementary charge
const double SI::mu0  = m_pi*4.0e-7;                // N A-2 vacuum permeability
const double SI::eps0 = 1.0/(SI::mu0*SI::c*SI::c);  // F m-1 vacuum permitivity

// Atomic constants
const double SI::me   = 9.1093897e-31;   // kg electron mass
const double SI::mp   = 1.6726231e-27;   // kg proton mass
const double SI::mn   = 1.6749286e-27;   // kg neutron mass

// Thermodynamic constants
const double SI::Kb   = 1.380658e-23;   // J K-1          Boltzmann constant
const double SI::L    = 6.0221367e23;   // mol-1          Avogadro constant
const double SI::R    = 8.314510;       // J mol-1 K-1    Gas constant

// Conversion factors
const double SI::eV   = SI::e;          // J  electron volt
const double SI::amu  = 1.6605402e-27;  // kg atomic mass unit
const double SI::atm  = 101325;         // Pa standard athmoshere
const double SI::cal  = 4.1840;         // J  thermodynamic calorie

