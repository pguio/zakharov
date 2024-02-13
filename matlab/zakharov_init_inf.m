function solver=zakharov_init_inf(solver)
% function solver=zakharov_init_inf(solver)

%
% $Id: zakharov_init_inf.m,v 1.10 2011/03/26 09:20:42 patrick Exp $
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

L = solver.L;
N = solver.N;

% Grid points xj in configuration space
solver.xj = [0:N-1]' * L / N;
solver.xj = solver.xj - L / 2.0;

% M < N/3 
solver.M = floor(N/3/(N/32));
M = solver.M;

% Grid points k in fourier space
solver.k = 2.0 * pi * [-N/2:N/2-1]' / L;

% Time initialisation
solver.ctime=0;

% Array allocation
solver.Ej = zeros(size(solver.xj));
solver.nj = zeros(size(solver.xj));
solver.dnj = zeros(size(solver.xj));

for k=1:solver.nb

	Emax = solver.Emax(k);
	m = solver.m(k);

% Eq. 31
	v = m *2*2*pi/L;
	solver.v(k) = v;

	Emin = 0;
	solver.Emin(k)=Emin;

	q = 1;
	solver.q(k)=q;

% Calculation of n0
% Eq. 27
  x = solver.xj;
	w = Emax/sqrt(2*(1-v^2))*x;

% Eq. 26
	F = Emax .* sech(w);
% Eq. 22
	n = -abs(F).^2/(1-v^2);
	n0 = -integrate(x(:),n(:))/L;
	solver.n0(k) = n0;

% Eq. 29
	u = v/2 + 2*n0/v - (Emax^2+Emin^2)/(v*(1-v^2));
	solver.u(k) = u;

	[E,n,dn]=analytic_inf_one_soliton(solver,k);


  solver.Ej = solver.Ej + E;
  solver.nj = solver.nj + n;
  solver.dnj = solver.dnj + dn;

end

solver.nk = FFT(solver.nj, N, M);
solver.nk(find(solver.k==0))=0;
solver.dnk = FFT(solver.dnj, N, M);
solver.dnk(find(solver.k==0))=0;
solver.Ek = FFT(solver.Ej, N, M);

if ~isfield(solver,'D'),
  solver.D = L;
end

if ~isfield(solver,'maxiter')
	h = solver.h;
	m = max(abs(solver.m(:)));
	v = max(abs(solver.v(:)));
	D = solver.D;
	solver.maxiter = round(D/v/h);
	solver.h = D/abs(v)/solver.maxiter;
	solver.maxtime = solver.maxiter *solver.h;
end

% Calculation of the cos and sin terms in Eq. 10 
k = solver.k;
nui = solver.nui;
h = solver.h;
ct = cos(sqrt(k.^2-nui.^2)*h);
st = repmat(h,size(ct));
ii = find(k~=0);
st(ii) = sin(sqrt(k(ii).^2-nui.^2)*h)./sqrt(k(ii).^2-nui.^2);
solver.ct = ct;
solver.st = st;

