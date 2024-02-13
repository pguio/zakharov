function solver=egauss_init(solver)
% function solver=egauss_init(solver)

%
% $Id: egauss_init.m,v 1.12 2011/03/26 09:20:41 patrick Exp $
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

% Gaussian wave packet
E0=sqrt(2*solver.W/solver.Delta/sqrt(2*pi)/(2*pi));
fprintf(1,'E0=%f\n', E0);
solver.Ek = E0*exp(-((solver.k-solver.k0)/solver.Delta).^2).* ...
	exp(i*2*pi*rand(size(solver.k)));
%E0=solver.W/(2*pi/solver.L)/integrate(solver.k,abs(solver.Ek).^2);
%solver.Ek=sqrt(E0)*solver.Ek;

% Random phase noise
solver.Ek = solver.Ek+solver.Erms*exp(i*2*pi*rand(size(solver.k)));
solver.Ek(find(solver.k==0)) = solver.E0;

% Random phase noise
%solver.nk = solver.nrms*exp(i*2*pi*rand(size(solver.k)));
solver.nk = solver.nrms*exp(i*2*pi*rand(size(solver.k)));
for i=1:solver.M,
	ip=find((solver.signed_index)==i);
	im=find((solver.signed_index)==-i);
	solver.nk(ip)=conj(solver.nk(im));
end
% Thermal fluctuations n(k=0)=0
solver.nk(find(solver.k==0)) = 0.0;

solver.dnk = zeros(size(solver.k));

% Set to zero k-modes larger than M
solver.Ek(solver.ikeq0)=0.0;
solver.nk(solver.ikeq0)=0.0;

% Normalise Fourier transform
solver.Ek=solver.Ek*solver.M;
solver.nk=solver.nk*solver.M;
solver.dnk=solver.dnk*solver.M;

% Fourier transform 
solver.Ej = IFFT(solver.Ek, solver.N, solver.M);
solver.nj = IFFT(solver.nk, solver.N, solver.M);
solver.dnj = IFFT(solver.dnk, solver.N, solver.M);

if 1
	subplot(121), plot(solver.k, abs(solver.nk).^2)
	set(gca,'xlim',[-20 20],'ylim',[0 5.5e-3]);
	subplot(111), plot(solver.k, abs(solver.Ek/solver.M).^2)
	set(gca,'xlim',[-20 20],'ylim',[0 5.5e-3]);
	%fprintf(1,'W=%f\n', (2*pi/solver.L)*integrate(solver.k,abs(solver.Ek).^2));
	fprintf(1,'W=%f\n', (2*pi)*integrate(solver.k,abs(solver.Ek/solver.M).^2));
	pause
end
