function solver=maxwell_param_init(solver)
% function solver=maxwell_param_init(solver)

%
% $Id: maxwell_param_init.m,v 1.9 2011/03/26 09:20:41 patrick Exp $
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

si=SIconstants;

% Ion mass
solver.mi=solver.mui*si.amu; 

solver.k = 2 * pi * [-solver.N/2:solver.N/2-1]' / solver.L;

solver.eta=(solver.Te+3*solver.Ti)/solver.Te;

% physical values
solver.wpe=sqrt(solver.ne*si.e^2/(si.eps0*si.me)); % [rad s^-1]
solver.ve=sqrt(si.kb*solver.Te/si.me); % [m s^-1]
solver.ve2=si.kb*solver.Te/si.me; % [m s^-1]^2
solver.cs=sqrt(si.kb*solver.Te/solver.mi); % [m s^-1]
solver.lambdae=sqrt(si.eps0*si.kb*solver.Te/(solver.ne*si.e^2)); % [m]
solver.lambdai=sqrt(si.eps0*si.kb*solver.Ti/(solver.ne*si.e^2)); % [m]

% normalisation constants
solver.tau=3/(2*solver.eta)*(solver.mi/si.me)/solver.wpe; % [s]

solver.chi=3/2*sqrt((solver.mi/si.me)/solver.eta)*solver.lambdae; % [m]

%solver.epsilon=4*sqrt(solver.eta*si.me/(3*solver.mi*si.eps0))* ...
%	sqrt(solver.ne*si.kb*solver.Te); % [V m^-1]
solver.epsilon=solver.eta*sqrt(si.me/solver.mi)* ...
	sqrt(4/3*solver.ne*si.kb*solver.Te/si.eps0); % [V m^-1]

solver.nu=4/3*solver.eta*(si.me/solver.mi); % [m^-3]


