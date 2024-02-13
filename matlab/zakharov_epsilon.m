function [h,epsilon]=zakharov_epsilon(solver)
% function [h,epsilon]=zakharov_epsilon(solver)
%
% Calculate epsilon given by Eq.33 of payne (1983)

%
% $Id: zakharov_epsilon.m,v 1.6 2011/03/26 09:20:42 patrick Exp $
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

if solver.nb~=1,
	error('zakharov_epsilon can handle only one soliton')
end

h = solver.h;
Eh = solver.Ej;

[E,n,dn] = analytic_one_soliton(solver,1);

if 0
	subplot(411), 
	plot(solver.xj,[real(E),real(Eh)]) 
	subplot(412), 
	plot(solver.xj,[imag(E),imag(Eh)]) 
	subplot(413), 
	plot(solver.xj,[abs(E),abs(Eh)]) 
	subplot(414),
	plot(solver.xj,[real(E-Eh),imag(E-Eh),abs(E-Eh)])
	legend('Re','Im','Abs')
	keyboard
end

X = solver.xj;

% Build the term of Eq. 33
dE = abs(E-Eh).^2;
E2 = abs(E).^2;

% Add the last point for the symmetry
X = [X(:); -solver.xj(1)];
dE = [dE(:); dE(1)];
E2 = [E2(:); E2(1)];

% Calculate the integrals and estimate epsilon
epsilon = sqrt(integrate(X, dE)/integrate(X, E2));

