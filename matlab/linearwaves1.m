function S = linearwaves(s,ispec,uspec,dspec)
% function S = linearwaves(s)

%
% $Id: linearwaves1.m,v 1.11 2011/03/26 09:20:41 patrick Exp $
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

if nargin == 0 | (nargin == 1 & ~isfield(s,'ito')),
	s.ito = 0;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'nk_stint')),
	s.nk_stint = 2;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'Ek_stint')),
	s.Ek_stint = 15;
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

if 0, plot_fields(s);
else plot_current_kspectra(s,kspec);
end
pause(1)


for j = 0:s.maxiter,
	s.j = j;
	[snk, sEk, svk] = source_update(s);

	nkh = s.nk.*s.soundDecay1+s.dnk.*s.soundDecay2+snk;
	s.dnk = s.dnk.*s.soundDecay4-(s.h/2)*s.k.^2.*(s.nk.*s.soundDecay4+nkh)+svk;
	s.nk = nkh;
	s.Ek = s.Ek.*s.langmuirDecay+sEk;

	kspec = kspectra_update(kspec, s);
	ispec = spectra_update(ispec, s);
	uspec = spectra_update(uspec, s);
	dspec = spectra_update(dspec, s);

	if mod(s.j,100) == 0 | s.j == s.maxiter,
		fprintf(1,'Iter = %6d/%6d%s', s.j, s.maxiter, char(8*ones(1,20)));

		if 0, plot_fields(s);
		else, plot_current_kspectra(s,kspec); end
		drawnow
	end
end
fprintf(1,'\n');

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
s.ikc = find([-s.N/2:s.N/2-1]' <= s.M & [-s.N/2:s.N/2-1]' >= -s.M)';

s.h = s.H/s.tau;
s.k = s.K*s.chi;

s.k0 = find(s.k == 0);

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
s.soundDecay3 = exp(-s.nui*s.h);
s.soundDecay4 = exp(-2*s.nui*s.h);
s.langmuirDecay = exp(-(i*s.k.^2+s.nue)*s.h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = source_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
si=SIconstants;
a2 = (s.K*s.lambdae).^2;
s.snk = sqrt(s.ne^2/s.Ne)*sqrt(1./(1+a2)./(1+a2+s.Te/s.Ti));
s.snk = s.snk/s.nu*s.N;
s.snk([s.pad0;s.k0]) = 0.0;
s.Dn = sqrt(2*s.nui.*s.k.^2./(s.k.^2+4*s.nui.^2)).*s.snk;
%subplot(111), plot(s.k, [s.snk, s.Dn]); pause

s.sEk = sqrt(s.ne^2/s.Ne)*(si.e/si.eps0)*s.lambdae./sqrt(2*(1+a2));
s.sEk = s.sEk/s.epsilon*s.N;
s.sEk([s.pad0;s.k0]) = 0.0;
s.DE = sqrt(s.nue).*s.sEk;
%subplot(111), plot(s.k, [s.sEk, s.DE]); pause

rand('twister', s.seed);

s.ip = find(s.k > 0 & [-s.N/2:s.N/2-1]' <= s.M);
s.im = find(s.k < 0 & [-s.N/2:s.N/2-1]' >= -s.M);

k = s.k(s.ikc);
nbk = length(k);
nui = s.nui(s.ikc);
signk = s.Dn(s.ikc);

h = [0:s.h/s.nk_stint:s.h];
dh = diff(h);
if s.ito, % Ito integral
	s.signks = repmat(sqrt(dh)',1,nbk).*repmat(signk',s.nk_stint,1).* ...
		exp(repmat(nui.',s.nk_stint,1).*repmat(h(1:end-1)',1,nbk));
else, % Stratonovich integral
	s.signks = repmat(sqrt(dh)',1,nbk).*repmat(signk',s.nk_stint,1).* ...
		exp(repmat(nui.',s.nk_stint,1).*repmat(0.5*(h(1:end-1)+h(2:end))',1,nbk));
end

if s.ito, % Ito integral
	s.sigvks = repmat(sqrt(dh)',1,nbk).*repmat(signk',s.nk_stint,1).* ...
		exp(repmat(2*nui.',s.nk_stint,1).*repmat(h(1:end-1)',1,nbk));
else, % Stratonovich integral
	s.sigvks = repmat(sqrt(dh)',1,nbk).*repmat(signk',s.nk_stint,1).* ...
		exp(repmat(2*nui.',s.nk_stint,1).*repmat(0.5*(h(1:end-1)+h(2:end))',1,nbk));
end

nue = s.nue(s.ikc);
sigEk = s.DE(s.ikc);

h = [0:s.h/s.Ek_stint:s.h];
dh = diff(h);
if s.ito, % Ito integral
	s.sigEks = repmat(sqrt(dh)',1,nbk).*repmat(sigEk',s.Ek_stint,1).* ...
		exp(repmat((i*k.^2+nue).',s.Ek_stint,1).*repmat(h(1:end-1)',1,nbk));
else, % Stratonovich integral
	s.sigEks = repmat(sqrt(dh)',1,nbk).*repmat(sigEk',s.Ek_stint,1).* ...
	exp(repmat((i*k.^2+nue).',s.Ek_stint,1).* ...
	repmat(0.5*(h(1:end-1)+h(2:end))',1,nbk));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [snk, sEk, svk] = source_update(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
snk = zeros(size(s.k));
u1 = rand(size(s.signks));
u2 = rand(size(s.signks));
wnk = sqrt(-2*log(u1)).*exp(i*2*pi*u2);
wnk(:,fix(end/2)+2:end) = fliplr(conj(wnk(:,1:fix(end/2))));
snk(s.ikc) = sum(s.signks.*wnk).';%.*s.soundDecay3(s.ikc);

sEk = zeros(size(s.k));
u1 = rand(size(s.sigEks));
u2 = rand(size(s.sigEks));
wEk = sqrt(-2*log(u1)).*exp(i*2*pi*u2);
sEk(s.ikc) = sum(s.sigEks.*wEk).';%.*s.langmuirDecay(s.ikc);

if 1
	svk = zeros(size(s.k));
else
	if 1
		u1 = rand(size(s.im));
		u2 = rand(size(s.im));
		wvk = zeros(size(s.k));
		wvk(s.im) = sqrt(s.h)*sqrt(-2*log(u2)).*exp(i*2*pi*u1);
		wvk(s.ip) = flipud(conj(wvk(s.im)));
		svk = s.Dn.*wvk;
	else
		svk = zeros(size(s.k));
		u1 = rand(size(s.signks));
		u2 = rand(size(s.signks));
		wvk = sqrt(-2*log(u2)).*exp(i*2*pi*u1);
		wvk(:,fix(end/2)+2:end) = fliplr(conj(wvk(:,1:fix(end/2))));
		svk(s.ikc) = sum(s.signks.*wvk).'.*s.soundDecay4(s.ikc);
	end
end

if 0,
	subplot(211), plot(s.k(s.ikc), abs(snk(s.ikc)));
	subplot(212), plot(s.k(s.ikc), abs(sEk(s.ikc)));
	pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = field_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%
s.nk = zeros(size(s.k));
s.dnk = zeros(size(s.k));
s.Ek = zeros(size(s.k));
if 1,
	s.nk = s.snk;

	s.dnk(s.im) = s.Dn(s.im).*s.k(s.im)./sqrt(2*s.nui(s.im)).* ...
		exp(i*acos( 4*s.nui(s.im).^2./(s.k(s.im).^2+4*s.nui(s.im).^2)));
	s.dnk(s.ip) = s.Dn(s.ip).*s.k(s.ip)./sqrt(2*s.nui(s.ip)).* ...
		exp(i*acos(-4*s.nui(s.ip).^2./(s.k(s.ip).^2+4*s.nui(s.ip).^2)));

	s.Ek = s.sEk;
else
	[s.nk, s.Ek, s.dnk] = source_update(s);
end

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
plot(s.k, s.nk.*conj(s.nk)*(s.nu/s.N)^2);
subplot(212), 
plot(s.k, s.Ek.*conj(s.Ek)*(s.epsilon/s.N)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_current_kspectra(s,spec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(211),
h=plot(spec.k(s.ikc), spec.nk2(s.ikc)/spec.nspec*(s.nu/s.N)^2, ...
	spec.k(s.ikc),s.snk(s.ikc).^2*(s.nu/s.N)^2);
set(h(2),'linewidth',2)
title('<|n_k(t)|^2>_t')
set(gca,'xticklabel',[]);
subplot(212),
h=plot(spec.k(s.ikc), spec.Ek2(s.ikc)/spec.nspec*(s.epsilon/s.N)^2, ...
	spec.k(s.ikc),s.sEk(s.ikc).^2*(s.epsilon/s.N)^2);
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


