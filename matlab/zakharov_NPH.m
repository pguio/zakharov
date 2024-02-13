function [h,N,P,H]=zakharov_NPH(solver)
% function [h,N,P,H]=zakharov_NPH(solver)

%
% $Id: zakharov_NPH.m,v 1.6 2011/03/26 09:20:42 patrick Exp $
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

x = solver.xj;
E = solver.Ej;
n = solver.nj;
dn = solver.dnj;

k = solver.k;
dnk = solver.dnk;
Ek = solver.Ek;

h = solver.h;

X = [x(:); -x(1)];

% Calculation of N (Eq. 34)

Ns = abs(E).^2;
Ns = [Ns(:); Ns(1)];

% integrate N
N = integrate(X, Ns);

% calculation of V (Eq. 37)
ii = find(k~=0);
VK = zeros(size(n));
VK(ii) = -dnk(ii)./(i*k(ii));
V =  real(IFFT(VK,solver.N,solver.M));

if 0
	dV = derivate(x, V);
	dV1 = real(IFFT(i*k.*FFT(V,solver.N,solver.M),solver.N,solver.M));
	plot(x, dn, x, -dV, x, -dV1)
	keyboard
end

% Calculation of P (Eq. 35)
%dE1 = derivate(x,real(E))+i*derivate(x,imag(E));
dE = IFFT(i*k.*Ek,solver.N,solver.M);

if 0
	subplot(311), plot(x,real(E),x,imag(E))
	subplot(312), plot(x,real(dE),x,imag(dE))
	subplot(313), plot(x,real(dE1),x,imag(dE1))
	fprintf(1,'mean diff= %e %e\n', ...
		mean(abs(real(dE-dE1))),mean(abs(imag(dE-dE1))))
	keyboard
end

Ps = i / 2.0 *(E.*conj(dE) - conj(E).*dE) + n.*V;
Ps = [Ps(:); Ps(1)];

% integrate P
P =  integrate(X, Ps);

% Calculation of H (Eq. 36)

Hs = abs(dE).^2 + n.* abs(E).^2 + 0.5* n.^2+ 0.5.*V.^2;
Hs = [Hs(:); Hs(1)];

% integrate H
H =  integrate(X, Hs);

