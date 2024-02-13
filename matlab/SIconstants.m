function si=SIconstants
% function si=SIconstants

%
% $Id: SIconstants.m,v 1.6 2011/03/26 09:20:41 patrick Exp $
%
% Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2.  of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%
% Universal constants
%%%%%%%%%%%%%%%%%%%%%
si.c=299792458.0;              % m s-1       speed of light
si.c_tol=0.0;

si.G=6.67259e-11;              % m3 /kg /s2  gravitation constant
si.G_tol=128.0e-6;

si.h=6.6260755e-34;            % J s         Planck constant
si.h_tol=0.60e-6;

si.hbar=si.h/(2*pi);           % J s         Planck constant/(2pi)
si.hba_tol=si.h_tol;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electromagnetic constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%
si.e=1.60217733e-19;           % C           elementary charge
si.e_tol=0.30e-6;

si.mu0=pi*4.0e-7;              % N /A2       vacuum permeability
si.mu0_tol=0.0;

si.eps0=1.0/(si.mu0*si.c^2);   % F /m        vacuum permitivity
si.epso_tol=0.0;

%%%%%%%%%%%%%%%%%%
% Atomic constants
%%%%%%%%%%%%%%%%%%
si.me=9.1093897e-31;           % kg          electron mass
si.me_tol=0.59e-6;

si.mp=1.6726231e-27;           % kg          proton mass
si.mp_tol=0.59e-6;

si.mn=1.6749286e-27;           % kg          neutron mass
si.mn_tol=0.59e-6;

si.alpha=7.29735308e-3;        % 1           fine-structure constant
si.alpha_tol=0.045e-6;

si.Ry=10973731.534;            % /m          Rydberg constant
si.Ry_tol=0.0012e-6;

si.a0=0.529177249e-10;         % m           Bohr radius
si.a0_tol=0.045e-6;

si.ge=2.002319304386;          % 1           electron g-factor
si.ge_tol=1.0e-11;

si.re=2.81794092e-15;          % m           classical electron radius
si.re_tol=0.13e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Thermodynamic constants
%%%%%%%%%%%%%%%%%%%%%%%%%
si.kb=1.380658e-23;            % J /K        Boltzmann constant
si.kb_tol=8.5e-6;

si.L=6.0221367e23;             % /mol        Avogadro constant
si.L_tol=0.59e-6;

si.R=8.314510;                 % J /mol /K   Gas constant
si.R_tol=8.4e-6;

%%%%%%%%%%%%%%%%%%%%
% Conversion factors
%%%%%%%%%%%%%%%%%%%%
si.eV=si.e;                    % J           electron volt
si.eV_tol=si.e_tol;

si.amu=1.6605402e-27;          % kg          atomic mass unit
si.amu_tol=0.59e-6;

si.atm=101325;                 % Pa          standard athmoshere
si.atm_tol=0.0;

si.cal=4.1840;                 % J           thermodynamic calorie
si.atm_tol=0.0;

%%%%%%%%%%%%%%%%%%%%%%%%
% Astronomical constants
%%%%%%%%%%%%%%%%%%%%%%%%
si.Ms=1.99e30;                 % kg          solar mass
si.Rs=6.96e8;                  % m           solar radius
si.Tcs=13.6e6;                 % K           solar central temperature
si.Tss=5760;                   % K           solar surface temperature
si.au=1.496e11;                % m           astronomical unit
si.pc=3.09e16;                 % m           parsec

