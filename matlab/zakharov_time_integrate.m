function solver=time_integrate_zakharov(solver)
% function solver=time_integrate_zakharov(solver)

%
% $Id: zakharov_time_integrate.m,v 1.21 2011/03/26 09:20:42 patrick Exp $
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

% from Payne et al., Numerical solution of the Zakharov equations,
% J. Comp. Phys., 50, 482-498, 1983

nui = solver.nui;
nue = solver.nue;
h = solver.h;
k = solver.k;

N = solver.N;
M = solver.M;

nk = solver.nk;
dnk = solver.dnk;
Ek = solver.Ek;

ct = solver.ct;
st = solver.st;

% Compute thermal source
solver=zakharov_sources(solver);

% Term (|E(0)|^2)_k in Eq. 10
%E=solver.edrive+solver.Ej;
Ek2 = FFT(solver.Ej.*conj(solver.Ej), N, M);

% Eq. 10
% Integrate nk from t to t+h
nkh = nk.*exp(-nui*h).*(ct+nui.*st) ...
	+ dnk.*exp(-nui*h).*st-(h/2)*k.^2.*Ek2.*exp(-nui*h).*st+solver.snk;

% Fourier tranform nk into nj needed for Eq. 15
njh = IFFT(nkh, N, M);

% Eq. 12
Eestk = Ek.*exp(-(i*k.^2+nue)*h);
% Fourier tranform Ek into Ej
Eestj = IFFT(Eestk, N, M);

% Term (n(0)E(0))_k in Eq. 13
%E=solver.edrive+solver.Ej;
nEk = FFT(solver.nj.*solver.Ej, N, M);

if 0
subplot(311), plot(k, abs(Ek2)), set(gca,'xlim',0.7*[min(k) max(k)]);
subplot(312), plot(k, abs(nkh)), set(gca,'xlim',0.7*[min(k) max(k)]);
subplot(313), plot(k, abs(nEk)), set(gca,'xlim',0.7*[min(k) max(k)]);
drawnow
pause
end

% Eq. 13
%nEestk = nEk .* exp(-(i*k.^2 -solver.domega + nue) * h);
nEestk = nEk.*exp(-(i*k.^2+nue)*h);
% Fourier tranform nEk into nEj
nEestj = IFFT(nEestk,N,M);

% Fourier transform of thermal source 
if length(solver.sEk)==N,
	sEj = IFFT(solver.sEk,N,M);
else
	sEj = 0.0;
end

% Eq. 15
Ejh = (Eestj-i*(h/2)*(nEestj+sEj))./(1.0+i*(h/2)*njh);
% Fourier tranform Ejh into Ekh
Ekh = FFT(Ejh,N,M);

% Term (n(k)E(k))_k in Eq. 19
Ekh2 = FFT(Ejh.*conj(Ejh), N, M);

% Eq. 19
dnkh = dnk.*exp(-2*nui*h)-(h/2)*k.^2.*((nk+Ek2).* ...
	exp(-2*nui*h)+nkh+Ekh2)+solver.sdnk;

solver.nk = nkh;
solver.dnk = dnkh;
solver.Ek = Ekh;

solver.nj = IFFT(nkh, N, M);
solver.dnj = IFFT(dnkh, N, M);
solver.Ej = IFFT(Ekh, N, M);

% If electron velocity distribution function defined
% diffusion to be considered
if isfield(solver,'Fej'),
	solver=diffusion_update(solver);
end

solver.ctime = solver.ctime+h;


