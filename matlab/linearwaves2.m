function S = linearwaves(s,ispec,uspec,dspec)
% function S = linearwaves(s)

%
% $Id: linearwaves2.m,v 1.3 2011/03/26 09:20:41 patrick Exp $
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


close all

if nargin == 0 | (nargin == 1 & ~isfield(s,'c')),
	s.c = 1;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'L')),
	s.L = 100;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'N')),
	s.N = 4096;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'H')),
	s.H = 1e-6/s.c;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'maxiter')),
	s.maxiter = 20480*s.c;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'ne')),
	s.ne = 5.0e11;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'Te')),
	s.Te = 3000;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'nuec')),
	s.nuec = 100;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'miamu')),
	s.miamu = 16;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'Ti')),
	s.Ti = 1000;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'nuic')),
	s.nuic = 0.1;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'seeds')),
	s.seeds = 1:10;
end

ispec.type='i';
if nargin == 0 | (nargin == 1 & ~isfield(ispec,'stride')),
	ispec.stride = 20*s.c;
end

if nargin == 0 | (nargin == 1 & ~isfield(ispec,'ns')),
	ispec.ns = 128;
end

if nargin == 0 | (nargin == 1 & ~isfield(ispec,'ik')),
	ispec.ik = 2198:5:2698;
end

uspec.type='u';
if nargin == 0 | (nargin == 1 & ~isfield(uspec,'stride')),
	uspec.stride = 1*s.c;
end

if nargin == 0 | (nargin == 1 & ~isfield(uspec,'ns')),
	uspec.ns = 128;
end

if nargin == 0 | (nargin == 1 & ~isfield(uspec,'ik')),
	uspec.ik = 2198:5:2698;
end

dspec.type='d';

if nargin == 0 | (nargin == 1 & ~isfield(dspec,'stride')),
	dspec.stride = 1*s.c;
end

if nargin == 0 | (nargin == 1 & ~isfield(dspec,'ns')),
	dspec.ns = 128;
end

if nargin == 0 | (nargin == 1 & ~isfield(dspec,'ik')),
	dspec.ik = 1400:5:1900;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for seed = s.seeds,
	s.seed = seed;
	fprintf(1,'seed=%d\n',s.seed);
	s = solve(s,ispec,uspec,dspec);
	if seed == s.seeds(1);
		S = s;
	else
		S.kspec.nk2 = S.kspec.nk2+s.kspec.nk2;
		S.kspec.Ek2 = S.kspec.Ek2+s.kspec.Ek2;
		S.ispec.skf2 = S.ispec.skf2+s.ispec.skf2;
		S.uspec.skf2 = S.uspec.skf2+s.uspec.skf2;
		S.dspec.skf2 = S.dspec.skf2+s.dspec.skf2;
	end
end
S.kspec.nk2 = S.kspec.nk2/length(s.seeds);
S.kspec.Ek2 = S.kspec.Ek2/length(s.seeds);
S.ispec.skf2 = S.ispec.skf2/length(s.seeds);
S.uspec.skf2 = S.uspec.skf2/length(s.seeds);
S.dspec.skf2 = S.dspec.skf2/length(s.seeds);

plot_spectra(S.dspec, S.ispec, S.uspec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=solve(s,ispec,uspec,dspec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = variable_init(s);
s = damping_init(s);
s = solver_init(s);
s = source_init(s);
s = field_init(s);

kspec = kspectra_init(s);
ispec = spectra_init(ispec,s); 
uspec = spectra_init(uspec,s); 
dspec = spectra_init(dspec,s); 


s.ikc = find([-s.N/2:s.N/2-1]' <= s.M & [-s.N/2:s.N/2-1]' >= -s.M)';
A = zeros(length(s.ikc),3,3);
for ik=s.ikc,
	k2 = s.k(ik)^2;
	nue = s.nue(ik);
	nui = s.nui(ik);
	A(ik,:,:) = expm([-(i*k2+nue),0,0; 0,0,1; 0,-k2,-2*nui]*s.h);
end

if 0, plot_fields(s);
else plot_current_kspectra(s,kspec);
end
pause(1)


for j = 0:s.maxiter,
	s.j = j;
	[s, snk, sEk, sdnk] = source_update(s);

	for ik=s.ikc,
		Y = squeeze(A(ik,:,:))*[s.Ek(ik); s.nk(ik); s.dnk(ik)];
		s.Ek(ik) = Y(1);
		s.nk(ik) = Y(2);
		s.dnk(ik) = Y(3);
	end
	s.Ek = s.Ek + sEk;
	s.nk = s.nk + snk;
	s.dnk = s.dnk + sdnk;

	kspec = kspectra_update(kspec, s);
	ispec = spectra_update(ispec, s);
	uspec = spectra_update(uspec, s);
	dspec = spectra_update(dspec, s);

	if mod(s.j,100) == 0 | s.j == s.maxiter,
		fprintf(1,'Iter %6.d/%6.d\n',s.j, s.maxiter);
		if 0, plot_fields(s);
		else, plot_current_kspectra(s,kspec); end
		drawnow
	end
end

[kspec,s] = kspec_normalise(kspec,s);
ispec = spec_normalise(ispec, s);
uspec = spec_normalise(uspec, s);
dspec = spec_normalise(dspec, s);

plot_spectra(dspec, ispec, uspec);

s.kspec = kspec;
s.ispec = ispec;
s.uspec = uspec;
s.dspec = dspec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = variable_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
si=SIconstants;
s.me = si.me;
s.mi = s.miamu*si.amu;
s.Ne = s.ne*s.L^3;

s.we = sqrt(s.ne*si.e^2/(si.eps0*s.me));
s.ve = sqrt(si.kb*s.Te/s.me);
s.eta = (s.Te+3*s.Ti)/s.Te;
s.lambdae = sqrt(si.eps0*si.kb*s.Te/(s.ne*si.e^2));

s.tau = 3/2*(s.mi/(s.me*s.eta))/s.we;
s.chi = 3/2*sqrt(s.mi/(s.me*s.eta))*s.lambdae;
s.epsilon = 4/sqrt(3)*s.eta*sqrt(s.me/s.mi)*sqrt(s.ne*si.kb*s.Te/si.eps0);
s.nu = 4/3*s.eta*(s.me/s.mi)*s.ne;

s.K = 2*pi* [-s.N/2:s.N/2-1]'/s.L;
s.M = floor(s.N/3);
s.pad0 = find([-s.N/2:s.N/2-1]' > s.M | [-s.N/2:s.N/2-1]' < -s.M);
s.h = s.H/s.tau;
s.k = s.K*s.chi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = damping_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s.nui = zeros(size(s.k));
s.nui = s.nuic/2*s.tau + sqrt(pi/8)* ...
	(sqrt(s.me/s.mi)+(s.Te/s.Ti)^1.5*exp(-s.Te/(2*s.Ti)-3/2))*abs(s.k);

s.nue = zeros(size(s.k));
ik = find(s.k~=0);
s.nue(ik) = s.nuec/2*s.tau + sqrt(pi/8)*(3/2)^4*(s.mi/(s.eta*s.me))^2.5./ ...
	abs(s.k(ik)).^3.*(1+4/3*s.k(ik).^2*s.eta*s.me/s.mi).* ...
	exp(-9/8*s.mi/(s.eta*s.me)./s.k(ik).^2-3/2);


%%%%%%%%%%%%&&%%%%%%%%%%%%%
function s = solver_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
frq = s.k.^2-s.nui.^2;
s.ct = ones(size(s.k));
s.st = repmat(s.h, size(s.k));
ip = find(frq>0);
s.ct(ip) = cos(sqrt(frq(ip))*s.h);
s.st(ip) = sin(sqrt(frq(ip))*s.h)./sqrt(frq(ip));
im=find(frq<0);
s.ct(im) = cosh(sqrt(-frq(im))*s.h);
s.st(im) = sinh(sqrt(-frq(im))*s.h)./sqrt(-frq(im));

s.soundDecay1 = exp(-s.nui*s.h).*(s.ct+s.nui.*s.st);
s.soundDecay2 = exp(-s.nui*s.h).*s.st;
s.soundDecay4 = exp(-2*s.nui*s.h);
s.langmuirDecay = exp(-(i*s.k.^2+s.nue)*s.h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = source_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
si=SIconstants;
a2 = (s.K*s.lambdae).^2;
s.snk = sqrt(s.ne^2/s.Ne)*sqrt(1./(1+a2)./(1+a2+s.Te/s.Ti))/s.nu*s.N;
s.snk(s.pad0) = 0.0;
s.signk = sqrt(4*s.nui.*s.k.^2./(s.k.^2+4*s.nui.^2)).*s.snk;
%subplot(111), plot(s.k, [s.snk, s.signk]); pause

if 0
a2 = (s.K*s.lambdae).^2;
else
a2 = (s.lambdae/s.chi).^2*ones(size(s.k));
a2(1) = 0;
end
s.sEk = sqrt((s.ne*si.kb*s.Te/si.eps0)/s.Ne)*sqrt(a2/2)/s.epsilon*s.N;
s.sEk(s.pad0) = 0.0;
%ii=find(s.k~=0); s.sEk(ii) = s.sEk(ii)./abs(s.k(ii));
s.sigEk = sqrt(2*s.nue).*s.sEk;
%subplot(111), plot(s.k, [s.sEk, s.sigEk]); pause

rand('state', s.seed);

s.ip = find(s.k > 0 & [-s.N/2:s.N/2-1]' <= s.M);
s.im = find(s.k < 0 & [-s.N/2:s.N/2-1]' >= -s.M);
s.i0 = find(s.k == 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s,snk, sEk, sdnk] = source_update(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = rand(size(s.im));
x2 = rand(size(s.im));
wnk = zeros(size(s.k));
wnk(s.im) = sqrt(s.h)*sqrt(-2*log(x2)).*exp(i*2*pi*x1)/sqrt(2);
%wnk(s.im) = sqrt(s.h)*exp(i*2*pi*x1);
wnk(s.ip) = flipud(conj(wnk(s.im)));
snk = s.signk.*wnk;

x1 = rand(size(s.k));
x2 = rand(size(s.k));
wEk = sqrt(s.h)*sqrt(-2*log(x2)).*exp(i*2*pi*x1)/sqrt(2);
%wEk = sqrt(s.h)*exp(i*2*pi*x1);
if 1
	sEk = s.sigEk.*wEk;
else
	sEk = s.sigEk.*exp(-(i*s.k.^2+s.nue)*s.h/2).*(wEk-s.wEk);
	s.wEk = wEk;
end

sdnk = zeros(size(s.k));

%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = field_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%
s.nk = zeros(size(s.k));
s.dnk = zeros(size(s.k));
s.Ek = zeros(size(s.k));
if 1
	x1 = rand(size(s.k(s.im)));
	x2 = rand(size(s.k(s.im)));
	s.nk(s.im) = s.snk(s.im).*exp(i*2*pi*x1);
	s.nk(s.ip) = flipud(conj(s.nk(s.im)));
	s.nk = s.snk.*exp(i*2*pi*(rand(size(s.k))));

	x1 = rand(size(s.k(s.im)));
	x2 = rand(size(s.k(s.im)));
	s.Ek(s.im) = s.sEk(s.im).*exp(i*2*pi*x2);
	x1 = rand(size(s.k(s.ip)));
	x2 = rand(size(s.k(s.ip)));
	s.Ek(s.ip) = s.sEk(s.ip).*exp(i*2*pi*x2);
end
s.wEk = zeros(size(s.k));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spec = kspectra_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec.k = s.K;
if 0
	spec.nk2 = zeros(size(spec.k));
	spec.Ek2 = zeros(size(spec.k));
	spec.nspec = 0;
else
	spec.nk2 = s.nk.*conj(s.nk);
	spec.Ek2 = s.Ek.*conj(s.Ek);
	spec.nspec = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spec = kspectra_update(spec,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec.nk2 = spec.nk2+s.nk.*conj(s.nk);
spec.Ek2 = spec.Ek2+s.Ek.*conj(s.Ek);
spec.nspec = spec.nspec+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spec,s] = kspec_normalise(spec,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec.nk2 = spec.nk2/s.N^2/spec.nspec*s.nu^2;
spec.Ek2 = spec.Ek2/s.N^2/spec.nspec*s.epsilon^2;

s.snk = s.snk/s.N*s.nu;
s.sEk = s.sEk/s.N*s.epsilon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spec = spectra_init(spec,s) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec.k = s.K(spec.ik);
spec.ts = spec.stride*s.H;
spec.fs = 1/spec.ts;
spec.f = spec.fs*[-spec.ns/2:spec.ns/2-1]/spec.ns;
spec.skf2 = zeros(length(spec.k), length(spec.f));
spec.skt = zeros(length(spec.k), length(spec.f));
spec.nspec = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spec = spectra_update(spec,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rem(s.j, spec.stride) == 0,
	it = rem(fix(s.j/spec.stride),spec.ns)+1;
	if strcmp(spec.type,'i'),
		spec.skt(:, it) = s.nk(spec.ik);
	else,
		spec.skt(:, it) = s.Ek(spec.ik);
	end
	if it == spec.ns,
		if strcmp(spec.type,'d'),
			for k=1:length(spec.ik),
				skf = fft(spec.skt(k,:));
				spec.skf2(k,:) = spec.skf2(k,:)+skf.*conj(skf);
			end
		else
			for k=1:length(spec.ik),
				skf = ifft(spec.skt(k,:));
				spec.skf2(k,:) = spec.skf2(k,:)+skf.*conj(skf);
			end
		end
		spec.nspec = spec.nspec+1;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spec = spec_normalise(spec,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if spec.nspec ~= 0,
	spec.skf2 = spec.skf2/s.N^2/spec.ns^2/spec.nspec;
end
if strcmp(spec.type,'i'),
	spec.skf2 = spec.skf2*s.nu^2;
	for k=1:length(spec.ik)
		spec.skf2(k,:) = fftshift(spec.skf2(k,:));
	end
else
	spec.skf2 = spec.skf2*s.epsilon^2;
	for k=1:length(spec.ik)
		spec.skf2(k,:) = s.K(spec.ik(k)).^2*fftshift(spec.skf2(k,:));
	end
	si = SIconstants;
	spec.skf2 = spec.skf2*s.me/s.mi*(si.eps0/(si.kb*s.Te))^2;
end
[f,k] = meshgrid(spec.f, spec.k);
spec.sk2 = integrate(f', spec.skf2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_spectra(dspec,ispec,uspec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(231), 
imagesc(dspec.f, dspec.k, dspec.skf2); 
colorbar
subplot(232), 
imagesc(ispec.f, ispec.k, ispec.skf2);
axis xy
colorbar
subplot(233), 
imagesc(uspec.f, uspec.k, uspec.skf2);
axis xy
colorbar
subplot(223)
plot(ispec.k, ispec.sk2)
subplot(224)
plot(-dspec.k, dspec.sk2, uspec.k, uspec.sk2)

%%%%%%%%%%%%%%%%%%%%%%%
function plot_fields(s)
%%%%%%%%%%%%%%%%%%%%%%%
subplot(211), 
plot(s.k, s.nk.*conj(s.nk));
subplot(212), 
plot(s.k, s.Ek.*conj(s.Ek));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_current_kspectra(s,spec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(211),
h=plot(spec.k(s.ikc), spec.nk2(s.ikc)/spec.nspec,spec.k(s.ikc),s.snk(s.ikc).^2);
set(h(2),'linewidth',2)
title('<|n_k(t)|^2>_t')
set(gca,'xticklabel',[]);
subplot(212),
h=plot(spec.k(s.ikc), spec.Ek2(s.ikc)/spec.nspec,spec.k(s.ikc),s.sEk(s.ikc).^2);
set(h(2),'linewidth',2)
title('<|E_k(t)|^2>_t')
xlabel('k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_kspectra(spec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(211),
plot(spec.k(s.ikc), spec.nk2(s.ikc));
title('<|n_k(t)|^2>_t')
set(gca,'xticklabel',[]);
subplot(212),
plot(spec.k(s.ikc), spec.Ek2(s.ikc));
title('<|E_k(t)|^2>_t')
xlabel('k')


