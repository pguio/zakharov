function [E,n,dn]=analytic_one_soliton(solver,k)
% function [E,n,dn]=analytic_one_soliton(solver,k)

%
% $Id: analytic_one_soliton.m,v 1.6 2011/03/26 09:20:41 patrick Exp $
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
M = solver.M;

t = solver.ctime;

Emax = solver.Emax(k);
x0 = solver.x0(k);
n0 = solver.n0(k);
q = solver.q(k);
v = solver.v(k);
u = solver.u(k);

% Eq. 24
Phi = v/2;

% Eq. 27
x = solver.xj + x0 - v*t;
x=mod(x+solver.L/2, solver.L)-solver.L/2;

w = Emax/sqrt(2*(1-v^2))*x;

% dn calculation
[sn,cn,dn] = ellipj(w,q^2);
% Eq. 26
F = Emax * dn;

% Eq. 20
x = solver.xj + x0 - u*t;
E = -i * F .* exp(i*Phi*x);

% Eq. 22
n = -abs(F).^2/(1-v^2) + n0;

% d(dn(u,k^2))/du = -k^2 sn(u,k^2) cn(u,k^2)
dF = Emax * v * Emax/sqrt(2*(1-v^2)) * q^2 * sn .* cn;
dn = -2 * F .* dF / (1-v^2);
% numerically
%x = solver.xj;
%dn = -v * Emax/sqrt(2*(1-v^2)) * derivate(x, n);

if 0
	dna = -2 * F .* dF / (1-v^2);
	dnn = -v * Emax/sqrt(2*(1-v^2)) * derivate(x, n);
	subplot(111)
	plot(x, dna, x, dnn);
	keyboard
end


