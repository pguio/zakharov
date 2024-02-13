function s = nyisub(s)
% function s = nyisub(s)
%
% ion acoustic wavenumber damping operator from linear theory

%
% $Id: nyisub.m,v 1.3 2011/03/26 09:20:41 patrick Exp $
%
% Copyright (c) 2001-2011 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

p3 = sqrt(pi/8);

temp1 = p3*(sqrt(1/s.mass_ratio)+(s.t_ratio^1.5)*exp(-.5*s.t_ratio-1.5));

nyi = zeros(size(1:s.n));

m=1:s.n2-1;
nyi(m+1)     = temp1*s.pi2ndx*m + s.nyic;
nyi(s.n+1-m) = nyi(m+1);

nyi(1) = 0;
nyi(s.n2+1)  = temp1*s.pi2ndx*s.n2 + s.nyic;

s.nyi = nyi(:);

