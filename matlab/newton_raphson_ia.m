function newton_raphson_ia(param)
% function newton_raphson_ia(param)

%
% $Id: newton_raphson_ia.m,v 1.10 2011/03/26 09:20:41 patrick Exp $
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

if ~exist('param','var'),
	ne=1e11;
	Te=2000;
	Ti=1000;
	mi=16;
else
	ne=param{1};
	Te=param{2};
	Ti=param{3};
	mi=param{4};
end

ks = linspace(9,40,40);
%ks = linspace(-2*pi*2048/100,2*pi*2047/100,4096);
wrs = zeros(size(ks));
wis = zeros(size(ks));
wrsa = zeros(size(ks));
wisa = zeros(size(ks));
wrsf = zeros(size(ks));
wisf = zeros(size(ks));

i = 1;
for k=ks,
	[wrs(i) wis(i)] = numeric(k,ne,Te,Ti,mi);
	[wrsa(i) wisa(i)] = analytic(k,ne,Te,Ti,mi);
	[wrsf(i) wisf(i)] = fitted(k,ne,Te,Ti,mi);
	i = i+1;
end
	
subplot(211), plot(ks, wrsa, ks, wrs)
ylabel('\omega_r [s-1]');
title(sprintf('ne=%.2e Te=%d Ti=%d mi=%d',ne,Te,Ti,mi))
legend('analytic','numeric',2)

if Te==2000 & Ti==1000,
	subplot(212), plot(ks, wisa, ks, wis, ks, wisf)
	xlabel('k [m-1]');
	ylabel('\omega_i [s-1]');
	legend('analytic','numeric','fitted Te/Ti=2',3)
else
	subplot(212), plot(ks, wisa, ks, wis)
	xlabel('k [m-1]');
	ylabel('\omega_i [s-1]');
	legend('analytic','numeric',3);
end


function [wr,wi]=numeric(k,ne,Te,Ti,mi)

[wr, wi]=analytic(k,ne,Te,Ti,mi);
%fprintf(1,'k = %.4g m-1 w = (%.4g, %.4g) rad/s\n', k, wr, wi);
res = drel(wr+sqrt(-1)*wi,k,ne,Te,Ti,mi);
fprintf(1,'k = %.4g m-1 w = (%.4g, %.4g) rad/s W(w) = (%.4g, %.4g)\n', ...
	k, wr, wi, real(res), imag(res));

if k==0.0,
	wr = 0.0;
	wi = 0.0;
	return;
end

dz=1e-5;
epsilon=1e-10;
w = wr+sqrt(-1)*wi;
i = 0;
while 1,
	dw = drel(w,k,ne,Te,Ti,mi)/( ...
		(drel(w*(1+dz),k,ne,Te,Ti,mi)-drel(w*(1-dz),k,ne,Te,Ti,mi))/(2*dz*w));
	w = w - dw;
	i = i+1;
	if abs(dw)<epsilon,
		break;
	end
end
res = drel(w,k,ne,Te,Ti,mi);
fprintf(1,'iter = %d w = (%.4g, %.4g) rad/s W(w) = (%.4g, %.4g)\n', ...
	i, real(w), imag(w), real(res), imag(res));
wr = real(w);
wi = imag(w);

function [wr, wi]=analytic(k,ne,Te,Ti,mi)
si = SIconstants;
cs = sqrt(si.kb*(Te+3*Ti)/(si.amu*mi));
wr = k*cs;
wi = -sqrt(pi/8)* ...
	(sqrt(si.me/(mi*si.amu))+(Te/Ti)^1.5*exp(-Te/(2*Ti)-3/2))*abs(wr);


function [wr, wi]=fitted(k,ne,Te,Ti,mi)
si = SIconstants;
cs = sqrt(si.kb*(Te+3*Ti)/(si.amu*mi));
wr = k*cs;
wi = -17.16*k.^2-329.178*abs(k);


function z=drel(om,k,ne,Te,Ti,mi)
if 1
  si = SIconstants;
  li2= si.kb*si.eps0/si.e^2*Ti/ne;
  le2 = si.kb*si.eps0/si.e^2*Te/ne;
  vi = sqrt(si.kb*Ti/(si.amu*mi));
  ve = sqrt(si.kb*Te/(si.me));
else
  li2 = 4.7623e+03*Ti/ne;
	le2 = 4.7623e+03*Te/ne;
	vi = 91.1839*sqrt(Ti/mi);
	ve = 3.8931e+03*sqrt(Te);
end
k2 = k*k;

z = 1+1/(k2*li2)*w1(om/abs(k*vi))+1/(k2*le2)*w1(om/abs(k*ve));

