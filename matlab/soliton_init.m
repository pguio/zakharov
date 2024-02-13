function solver=soliton_init(solver)
% function solver=soliton_init(solver)

%
% $Id: soliton_init.m,v 1.5 2011/03/26 09:20:42 patrick Exp $
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

solver.Ej=zeros(size(solver.k));
solver.nj=zeros(size(solver.k));
solver.dnj=zeros(size(solver.k));

for k=1:solver.nb

	Emax = solver.Emax(k);
	m = solver.m(k);

% Eq. 31
	v = m * 2 * 2 * pi / solver.L;
	solver.v(k) = v;

% Emin Calculation
	if ~isfield(solver,'Emin'),
% Eq. 30
		q = inv_elliptke(solver.L*Emax/(2*sqrt(2*(1-v^2))));
		Emin = Emax*sqrt(1-q^2);
		solver.Emin(k) = Emin;
	else
% Eq. 28
		Emin = solver.Emin(k);
%		q = sqrt(Emax^2-Emin^2)/Emax;
		q = 1-0.5*(Emin/Emax)^2;
	end
	solver.q(k) = q;

% Calculation of n0
% Eq. 27
	x = solver.xj;
	w = Emax/sqrt(2*(1-v^2))*x;
% dn calculation
	[sn,cn,dn] = ellipj(w,q^2);
% Eq. 26
	F = Emax * dn;
% Eq. 22
	n = -abs(F).^2/(1-v^2);
	n0 = -integrate(x(:),n(:))/solver.L;
	solver.n0(k) = n0;

% Eq. 29
	u = v/2 + 2*n0/v - (Emax^2+Emin^2)/(v*(1-v^2));
	solver.u(k) = u;

	[E,n,dn]=analytic_one_soliton(solver,k);

	solver.Ej = solver.Ej + E;
	solver.nj = solver.nj + n;
	solver.dnj = solver.dnj + dn;

end

% Fourier transform 
solver.nk = FFT(solver.nj, solver.N, solver.M);
solver.nk(find(solver.k==0))=0;
solver.dnk = FFT(solver.dnj, solver.N, solver.M);
solver.dnk(find(solver.k==0))=0;
solver.Ek = FFT(solver.Ej, solver.N, solver.M);

if ~isfield(solver,'D'),
	solver.D = solver.L;
end

if ~isfield(solver,'maxiter')
	h = solver.h;
	v = max(abs(solver.v(:)));
	D = solver.D;
	solver.maxiter = round(D/(v*h));
	solver.h = D/abs(v)/solver.maxiter;
end
solver.maxtime = solver.maxiter *solver.h;

