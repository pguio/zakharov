function testdensfluct(seed,l0_L,lm_l0)
% function testdensfluct(seed, l0_L, lm_l0)
% 
% examples
% 
% testdensfluct(2,1,[1/2, 1/3/sqrt(2), 1/3/2])
% testdensfluct(2,1/2,1/4)

%
% $Id: testdensfluct.m,v 1.6 2011/03/26 09:20:42 patrick Exp $
%
% Copyright (c) 2010-2011 Patrick Guio <patrick.guio@gmail.com>
%
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

if ~exist('seed'),
  seed = 1;
end

close all

figure(1)

L0  = 90;
N0  = 8192;
dn0 = 2;

L  = [1, 2, 4]*L0;
N  = [1, 2, 4]*N0;

if ~exist('l0_L'),
  l0 = [1, 2, 4]*L0;
else
  l0 = l0_L.*L;
end

if ~exist('lm'),
  lm = [1/2, 1/3/sqrt(2), 1/3/2].*l0;
%  lm = 1/4*l0 ;
else
  lm = lm_l0.*l0;
end

for i = 1:length(L),

  figure(i)
  [xdnj{i},kdnk{i}] = init(L(i),N(i),dn0,l0(i),lm(i),seed);

end

figure

subplot(211)
plot(xdnj{1}{:},'-',xdnj{2}{:},'-',xdnj{3}{:},'-')
set(gca,'xlim',[-1 1]*max(L/2));

subplot(212)
plot(kdnk{1}{:},'-o',kdnk{2}{:},'-o',kdnk{3}{:},'-o')
set(gca,'xlim',[-1 1]*3*2*pi/min(lm))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xdnj,kdnk] = init(L,N,dn0,l0,lm,seed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'L = %6.2f N = %5d l0 = %6.2f lm = %6.2f\n', L, N, l0, lm);

% Grid points xj in configuration space
x = [0:N-1]' * L / N - L / 2.0;

% Grid points k in fourier space
k = 2 * pi * [-N/2:N/2-1]' / L;

in = find(k<0);
ip = find(k>0);

dnk = zeros(size(k));
dnk(k==0) = 0.0;

k0 = 2*pi/l0;
km = 2*pi/lm;
dk = km-k0;

dnk(in) = exp(-((-k(in)-k0)/dk).^2);
dnk(ip) = exp(-((k(ip)-k0)/dk).^2);

RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed));

phase = zeros(size(k));
phase(in) = i*2*pi*randn(size(k(in)));
phase(ip) = flipud(conj(phase(in(2:end))));

dnk = dnk.*exp(phase);

dnj = ifft(ifftshift(dnk));

fprintf(1,'imaginary part max %.4g\n', max(abs(imag(dnj)))); % = 0 0 is real exactly
if  max(abs(imag(dnj))) == 0,
dnj = real(dnj);
else
fprintf(1,'not exactly zero!\n');
end
dnj = dnj./max(abs(real(dnj)))*dn0;

dnk0 = dnk;
dnk = fftshift(fft(dnj));

subplot(211), plot(x,[real(dnj),imag(dnj)])
set(gca,'xlim',[min(x),max(x)])

subplot(212), plot(k,[abs(dnk),real(dnk),imag(dnk)],'-o')
set(gca,'xlim',[-1 1]*(3*km))



xdnj = {x, dnj};
kdnk = {k, abs(dnk)};

