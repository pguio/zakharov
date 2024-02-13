function x = ode1a(t,x0)
% function x = ode1a(t,x0)

%
% $Id: ode1a.m,v 1.3 2011/03/26 09:20:41 patrick Exp $
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

global nui k sign2 sigv2

if k~=0,
	A = [0,0,1;0,-4*nui,-k^2;-2*k^2,2,-2*nui];
	[v,d] = eig(A);
	r = [d(1,1),d(2,2),d(3,3)];

	l = 2*[ -nui+sqrt(nui^2-k^2), -nui, -nui-sqrt(nui^2-k^2) ];
	g1 = [1;k^2-2*sqrt(nui^2-k^2)*(nui-sqrt(nui^2-k^2));...
		-2*(nui-sqrt(nui^2-k^2))];
	g2 = [1;k^2;-2*nui];
	g3 = [1;k^2+2*sqrt(nui^2-k^2)*(nui+sqrt(nui^2-k^2));...
		-2*(nui+sqrt(nui^2-k^2))];
	g = [g1/norm(g1), g2/norm(g2), g3/norm(g3)];

	b = [sign2;sigv2;0];

	if 0
		l = r;
		g = v;
	end

	c = g\(x0+A\b);

	x = real(repmat(-(A\b).',length(t),1)+...
		c(1)*exp(l(1)*t)*g(:,1).'+...
		c(2)*exp(l(2)*t)*g(:,2).'+...
		c(3)*exp(l(3)*t)*g(:,3).');

else

	x = zeros(length(t),3);
	x(:,1) = x0(2)/(4*nui^2)*(exp(-4*nui*t)-1)+...
		sigv2/(4*nui^2)*(t+1/(4*nui)*(exp(-4*nui*t)-1))-...
		1/(2*nui)*(x0(3)+x0(2)/nui-sigv2/(2*nui^2))*(exp(-2*nui*t)-1)+...
		sign2*t+x0(1);
	x(:,2) = exp(-4*nui*t)*x0(2)+sigv2/(4*nui)*(1-exp(-4*nui*t));
	x(:,3) = -x0(2)/nui*exp(-4*nui*t)+...
		sigv2/(4*nui^2)*(1+exp(-4*nui*t))+(x0(3)+x0(2)/nui-...
		sigv2/(2*nui^2))*exp(-2*nui*t);
end

