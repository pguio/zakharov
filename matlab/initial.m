function [n0,dn0,e0]=initial(s)
% function [n0,dn0,e0]=initial(s)
% 
% initial profiles for n, dn/dt and E
% envelope soliton/random fluctuations

%
% $Id: initial.m,v 1.5 2011/03/26 09:20:41 patrick Exp $
%
% Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2. of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

if s.profile == 1, % soliton

  e0 = zeros(size([1:s.n]'));
  n0 = zeros(size([1:s.n]'));
  dn0 = zeros(size([1:s.n]'));

  v2 = s.v/2;
	pi2i = pi*s.im;
  w = sqrt(v2^2-s.omega);
	w2 = w^2;
	z = sqrt(2*(1-s.v^2));

	j = [1:s.n]';

	x = -12+s.dx*(j-1);
	zc1 = -s.im*v2*(x-s.x0);
	zc2 = s.im*v2*(x+s.x0);
	arg1 = w*(x-s.x0);
	arg2 = w*(x+s.x0);
	t1 = 1./cosh(arg1);
	t2 = 1./cosh(arg2);
	n0 = -2*w2*(t1.^2+t2.^2);
	dn0 = -4*s.v*w*w2*(-tanh(arg1).*t1.^2+tanh(arg2).*t2.^2);
	e0 = z*w*(exp(zc1).*t1+exp(zc2).*t2);

elseif s.profile == 2, % random fluctuations

  rand('twister',s.seed);

  fte0 = zeros(size([1:s.n]'));
  ftn0 = zeros(size([1:s.n]'));
  ftdn0 = zeros(size([1:s.n]'));

	j=[1:s.mmax]';
	
	fte0(j+1)      = s.fact1*exp(i*2*pi*rand(size(j)));
	fte0(s.n+1-j)  = s.fact1*exp(i*2*pi*rand(size(j)));
	ftn0(j+1)      = s.fact2*exp(i*2*pi*rand(size(j)));
	ftn0(s.n+1-j)  = conj(ftn0(j+1));
	ftdn0(j+1)     = 0;
	ftdn0(s.n+1-j) = 0;

	fte0(1) = 0; % no noise at k=0
	ftn0(1) = 0; % no noise at k=0

	ftn0(s.n2+1) = 2*real(ftn0(s.n2+1)); % nyquest mode is real 

  % zero padding of aliased modes when k=0 is first element
	fte0(s.zm)  = 0;
	ftn0(s.zm)  = 0;
	ftdn0(s.zm) = 0;

	e0  = ifft(fte0);
	n0  = ifft(ftn0);
	dn0 = ifft(ftdn0);


  if 0
	e02=fftshift(abs(fte0).^2);
	n02=fftshift(abs(ftn0).^2);
	e02(s.klim(1):s.klim(1)+1)
	e02(s.klim(2)-1:s.klim(2))
	subplot(211), plot(e02), set(gca,'xlim',s.klim)
	subplot(212), plot(n02), set(gca,'xlim',s.klim)
	pause
	end

end

