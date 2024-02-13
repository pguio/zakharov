function [out1,out2,out3] = ode1(t,y)
% function [out1,out2,out3] = ode1(t,y)

%
% $Id: ode1.m,v 1.3 2011/03/26 09:20:41 patrick Exp $
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

if length(t) == 0 			% return default tspan, y0 and options
  out1 = [0; 25];
  out2 = [0; 0; 0];
  out3 = [];
  return;
end

dy = zeros(size(y)); 			% preallocate vector dy

global nui k sign2 sigv2

dy(1) = y(3)+sign2;
dy(2) = -4*nui*y(2)-k^2*y(3)+sigv2;
dy(3) = -2*k^2*y(1)+2*y(2)-2*nui*y(3);

out1 = dy;
