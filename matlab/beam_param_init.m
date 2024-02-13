function solver=beam_param_init(solver)
% function solver=beam_param_init(solver)

%
% $Id: beam_param_init.m,v 1.8 2011/03/26 09:20:41 patrick Exp $
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

% Beam density
% Two possibilties : 
%     - fixed gamma0, nb_n0 is deduced
%     - fixed nb_n0, gamma0 is deduced
if isfield(solver,'gamma0'),
	solver.nb_n0=solver.gamma0/(sqrt(pi)*solver.tau*solver.wpe);
else,
	solver.gamma0=solver.nb_n0*sqrt(pi)*solver.tau*solver.wpe;
end

% Beam energy
% Two possibilties : 
%     - fixed k0 and vb is deduced
%     - fixed vb and k0 is deduced
if isfield(solver,'k0'),
% ici k0 est fixe et on calcul vb en fonction
  solver.vb=solver.wpe/solver.k0*solver.chi;
else
% Fixed vb 
  solver.k0=solver.wpe/solver.vb*solver.chi;
end

% Beam energy in eV
solver.Eb=0.5*si.me*solver.vb^2/si.eV;

