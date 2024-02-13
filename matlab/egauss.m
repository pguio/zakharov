function solver=egauss
% function solver=egauss

%
% $Id: egauss.m,v 1.12 2011/03/26 09:20:41 patrick Exp $
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

% spatial grid size
kmax=20;
solver.M=600; T=175; 
%solver.M=1200; T=300;
%solver.M=2400; T=550;
solver.N=pow2(nextpow2(solver.M*3));
% spatial period
solver.L=2*pi*solver.M/kmax;
% time increment
solver.h=1e-3;
% damping
solver.nui=0;
solver.nue=0;

solver.type='egauss';
solver.k0=8;
solver.W=0.04;
solver.Delta=1.0;
solver.E0=0.0; % initial E-field power
solver.Erms=1e-3; % initial E-field fluctuation
solver.Erms=1e-6; % initial E-field fluctuation
solver.nrms=1e-6; % initial density fluctuation
solver.nrms=1e-8; % initial density fluctuation

solver.display=1;
solver.mdisplay=100;
solver.displaytype='k';

solver.save=1;
solver.msave=1000;

solver.maxiter = T/solver.h;
solver.maxtime = solver.maxiter*solver.h;


