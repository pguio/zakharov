function solver=beam
% function solver=beam

%
% $Id: ekbeam.m,v 1.6 2011/03/26 09:20:41 patrick Exp $
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

% spatial period
%solver.L=15*pi;
solver.L=47.1238898038468932;
solver.L=80;
% spatial grid size
solver.N=4096;
% time increment
solver.h=1e-3;

% initial E-field power
solver.E0=0.0; 
% initial E-field fluctuation
solver.Erms=1e-4; 

solver.edrive=0.0; % electric pump field
solver.domega=0.0; % frequency mismatch

% damping
solver.lhanssen=1;
solver.nuic=0.01; % ion collisions [s^-1]
solver.nuec=2.4;  % electron collisions [s^-1]
solver.nui0=0.3;  % dimensionless ion Landau damping 

% temperatures
solver.Te=3000.0;
solver.Ti=1200.0;

% beam parameter
solver.gamma0=2.0*5;
%solver.k0=-22.0;
solver.k0=-10.0;
solver.dvb_vb=0.3;

% plasma parameters
solver.ne=7.e10;
solver.mui=16; % Molecular oxygen number

% si k0 est libre decomenter la ligne suivante 
if ~isfield(solver,'k0'),
	solver.vb=6.6235e6;
end

% Initialise background parameters
solver = maxwell_param_init(solver);
% Initialise beam parameters
solver = beam_param_init(solver);
% Collision frequencies
solver = collisions(solver);
% Initialise diffusion
%solver = diffusion_init(solver);

solver.type='tfluct';

solver.display=1;
solver.mdisplay=100;
solver.displaytype='k';

solver.save=1;
solver.msave=150;

solver.maxiter = 4500;
solver.maxtime = solver.maxiter*solver.h;

solver.av_start=1;
solver.av_end=solver.maxiter;

% Wave vectors to consider to calculate IS spectra
% Radar parameters
si=SIconstants;
f=224e6; 
solver.is_kmin=2*2*pi*f/si.c; % [m-1]
f=931e6;
solver.is_kmax=2*2*pi*f/si.c; % [m-1]
solver.is_nbk=12;

solver.is_start=1;
solver.is_end=solver.maxiter;
solver.is_stride=1;
%solver.is_nfft=1024;


