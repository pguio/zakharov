function solver=efluct0_1
% function solver=efluct0_1

%
% $Id: efluct0_1.m,v 1.7 2011/03/26 09:20:41 patrick Exp $
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
solver.L=160*pi;
% spatial grid size
solver.N=512;
% time increment
solver.h=1e-1;

solver.mass_ratio=40000;
solver.t_ratio=2;

solver.type='efluct';
solver.E0=1e-1; % initial E-field power
solver.Erms=1e-4; % initial E-field fluctuation

solver.edrive=0.0; % electric pump field
solver.domega=0.0; % frequency mismatch

k = 2 * pi * [-solver.N/2:solver.N/2-1]' / solver.L;

% damping
if 0
	solver.nuic=0;
	solver.nui=sqrt(pi/8)*(sqrt(1/solver.mass_ratio)+solver.t_ratio^1.5* ...
  	exp(-0.5*solver.t_ratio-1.5))*abs(k)+solver.nuic;

	solver.ecoll=0.5;
	eta=1+3/solver.t_ratio;
	s1=sqrt(pi/8)*1.5^4*exp(-1.5)*(solver.mass_ratio/eta)^2.5;
	s2=9.0*solver.mass_ratio/(8*eta);
	k2=k.*[-solver.N/2:solver.N/2-1]';
	solver.nue=zeros(size(k))+solver.ecoll/2;
	ii=find(k~=0);
	solver.nue(ii)=(s1./(k2(ii).*abs(k(ii)))).*exp(-s2./k2(ii))+solver.ecoll/2;
else
	solver.nui=0;
	solver.nue=0;
end


solver.display=1;
solver.mdisplay=100;
solver.displaytype='k';

solver.save=1;
solver.msave=fix(200000/1000);

solver.maxiter = 200000;
solver.maxtime = solver.maxiter*solver.h;

solver.av_start=20000;
solver.av_end=solver.maxiter;

solver.is_ks=solver.N/2+[10:20:150]+1;
solver.is_start=20000;
solver.is_start=1;
solver.is_end=solver.maxiter;
solver.is_stride=1;
%solver.is_nfft=1024;


