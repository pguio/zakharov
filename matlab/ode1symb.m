function ode1symb
% function ode1symb

%
% $Id: ode1symb.m,v 1.2 2011/03/26 09:20:41 patrick Exp $
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

s=dsolve('Dx=z+sn^2','Dy=-4*nui*y-k^2*z+sv^2','Dz=-2*k^2*x+2*y-2*nui*z',...
	'x(0)=x0','y(0)=y0','z(0)=z0');


s=dsolve('Dx=z+sn^2','Dy=-4*nui*y+sv^2','Dz=2*y-2*nui*z',...
	'x(0)=x0','y(0)=y0','z(0)=z0');

simplify(s.x)
