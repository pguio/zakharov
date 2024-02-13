function w=disproot(k,w,ne,Te,Ti,mi)
% function w=disproot(k,w,ne,Te,Ti,mi)

%
% $Id: disproot.m,v 1.3 2011/03/26 09:20:41 patrick Exp $
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

w = numeric(k,w,ne,Te,Ti,mi);

function w=numeric(k,w,ne,Te,Ti,mi)

res = drel(w,k,ne,Te,Ti,mi);
%fprintf(1,'k = %.4g m-1 w = (%.4g, %.4g) rad/s W(w) = (%.4g, %.4g)\n', ...
%	k, real(w), imag(w), real(res), imag(res));

if k==0.0,
	w = 0.0;
	return;
end

dz=1e-5;
epsilon=1e-10;
i = 0;
while 1,
	dw = drel(w,k,ne,Te,Ti,mi)/( ...
		(drel(w*(1+dz),k,ne,Te,Ti,mi)-drel(w*(1-dz),k,ne,Te,Ti,mi))/(2*dz*w));
	w = w - dw;
	i = i+1;
	if abs(dw)<epsilon | i>100,
		break;
	end
end
res = drel(w,k,ne,Te,Ti,mi);
%fprintf(1,'iter = %d w = (%.4g, %.4g) rad/s W(w) = (%.4g, %.4g)\n', ...
%	i, real(w), imag(w), real(res), imag(res));

function [wr, wi]=analytic(k,ne,Te,Ti,mi)
si = SIconstants;
cs = sqrt(si.kb*(Te+3*Ti)/(si.amu*mi));
wr = k*cs;
wi = -sqrt(pi/8)* ...
	(sqrt(si.me/(mi*si.amu))+(Te/Ti)^1.5*exp(-Te/(2*Ti)-3/2))*abs(wr);


function z=drel(om,k,ne,Te,Ti,mi)
si = SIconstants;
li2= si.kb*si.eps0/si.e^2*Ti/ne;
le2 = si.kb*si.eps0/si.e^2*Te/ne;
vi = sqrt(si.kb*Ti/(si.amu*mi));
ve = sqrt(si.kb*Te/(si.me));
k2 = k*k;

z = 1+1/(k2*li2)*w1(om/abs(k*vi))+1/(k2*le2)*w1(om/abs(k*ve));

