function c=inv_elliptke(x)
% function c=inv_elliptke(x)

%
% $Id: inv_elliptke.m,v 1.2 2011/03/26 09:20:41 patrick Exp $
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

a=0;
b=1-eps^2;
maxiter=200;

ya=ellipke(a^2)-x;
yb=ellipke(b^2)-x;

if ya*yb>0 
	error('inv_elliptke not well condtioned...')
end

i=0;
while abs(ya-yb)>1e-6 & i<maxiter,
	c=0.5*(a+b);
	yc=ellipke(c^2)-x;
	if ya*yc<0,
		b=c;
		yb=ellipke(b^2)-x;
	else
		a=c;
		ya=ellipke(a^2)-x;
	end
%	fprintf(1,'.');
	i=i+1;
end
c=0.5*(a+b);


