function ode2symb
% function ode2symb

%
% $Id: ode2symb.m,v 1.3 2011/03/26 09:20:41 patrick Exp $
%
% Copyright (c) 2001-2011 Patrick Guio <patrick.guio@gmail.com>
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

syms nue nui k t

A = sym([0,1;-k^2,-2*nui]*t)
S = expm(A);
pretty(simplify(S))

pause

A = sym([-nue-sqrt(-1)*k^2,0,0;0,0,1;0,-k^2,-2*nui]*t)
S = expm(A);
simplify(S)


A = sym([-nue-sqrt(-1)*k^2,0,0;0,0,1;0,-k^2,-2*nui]*t)


