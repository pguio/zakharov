function [E,n,dn]=analytic_inf_one_soliton(solver,k)
% function [E,n,dn]=analytic_inf_one_soliton(solver,k)

%
% $Id: analytic_inf_one_soliton.m,v 1.7 2011/03/26 09:20:41 patrick Exp $
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
q = solver.q;
v = solver.v;
u = solver.u;

% Eq. 24
Phi = v/2;

% Eq. 27
x = solver.xj + x0 - v*t;
x=mod(x+solver.L/2, solver.L)-solver.L/2;

w = Emax/sqrt(2*(1-v^2))*x;

% Eq. 26
F = Emax .* sech(w);

x = solver.xj + x0 - u*t;
E = -i * F .* exp(i*Phi*x);

% Eq. 22
n = -abs(F).^2/(1-v^2) + n0;

dF = - Emax * v * Emax/sqrt(2*(1-v^2)) .* sinh(w) .* sech(w).^2;
dn= 2 / (1 - v^2) * F .* dF;
%x = solver.xj;
%dn1 = -v * Emax/sqrt(2*(1-v^2)) * derivate(x, n);

if 0
	plot(x, dn, x, dn1);
end

