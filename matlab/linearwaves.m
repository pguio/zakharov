function s = linearwaves(s)
% function s = linearwaves(s)

%
% $Id: linearwaves.m,v 1.37 2011/03/26 09:20:41 patrick Exp $
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

if nargin == 0 | (nargin == 1 & ~isfield(s,'L_r')),
	s.L_r = 100;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'N')),
	s.N = 4096;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'h_r')),
	s.h_r = 1e-6/s.c; % time step [s]
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'maxiter')),
	s.maxiter = 5120*s.c;
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
	s.seeds = 1;
end

s.ispec.type='i';
if nargin == 0 | (nargin == 1 & ~isfield(s.ispec,'stride')),
	s.ispec.stride = 20*s.c;
end
if nargin == 0 | (nargin == 1 & ~isfield(s.ispec,'ns')),
	s.ispec.ns = 128;
end
if nargin == 0 | (nargin == 1 & ~isfield(s.ispec,'ik')),
	%s.ispec.ik = 2198:5:2698; % valid for N=4096
	s.ispec.ik = [fix(s.N/2)+1:fix(s.N/6/100):fix(s.N/2)+1+fix(s.N/6)];
end

s.uspec.type='u';
if nargin == 0 | (nargin == 1 & ~isfield(s.uspec,'stride')),
	s.uspec.stride = 1*s.c;
end
if nargin == 0 | (nargin == 1 & ~isfield(s.uspec,'ns')),
	s.uspec.ns = 128;
end
if nargin == 0 | (nargin == 1 & ~isfield(s.uspec,'ik')),
	%s.uspec.ik = 2198:5:2698; % valid for N=4096
	s.uspec.ik = [fix(s.N/2)+1:fix(s.N/6/100):fix(s.N/2)+1+fix(s.N/6)];
end

s.dspec.type='d';
if nargin == 0 | (nargin == 1 & ~isfield(s.dspec,'stride')),
	s.dspec.stride = 1*s.c;
end
if nargin == 0 | (nargin == 1 & ~isfield(s.dspec,'ns')),
	s.dspec.ns = 128;
end
if nargin == 0 | (nargin == 1 & ~isfield(s.dspec,'ik')),
	%s.dspec.ik = 1400:5:1900; % valid for N=4096
	s.dspec.ik = fliplr([fix(s.N/2)-1:-fix(s.N/6/100):fix(s.N/2)-1-fix(s.N/6)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s.ns = length(s.seeds);

s = variable_init(s);
s = damping_init(s);
s = solver_init(s);

fprintf(1,'Ion-acoustic and Langmuir linear wave equations\n\n');
fprintf(1,'\nDimensions\n');
fprintf(1,'maxiter=%d N=%d M=%d\n', s.maxiter, s.N, s.M);
fprintf(1,'L=%g [m] (%g)\nh=%g [s] (%g)\n', s.L_r, s.L_r/s.chi, s.h_r, s.h_r/s.tau);
fprintf(1,'\nPlasma parameters and derived\n')
fprintf(1,'ne=%g [m-3]\nTe=%g [K] Ti=%g [K]\nmi=%g [amu]\nnuic=%g [s-1] nuec=%g [s-1]\n', ...
	s.ne, s.Te, s.Ti, s.miamu, s.nuic, s.nuec);
fprintf(1,'fpe=%.3g [MHz] ve=%.3g [m s-1] lambdae=%.3g [m]\n', s.we/(2*pi)*1e-6, s.ve, s.lambdae);
fprintf(1,'\nNormalisation constants\n');
fprintf(1,'tau=%g [s]\nchi=%g [m]\nepsilon=%g [V m-1]\nnu=%g [m-3]\n', ...
	s.tau, s.chi, s.epsilon, s.nu);

for is = 1:s.ns,
	s.is = is;
	fprintf(1,'seed=%d/%d\n',s.is, s.ns);
	s = solve(s);
end

plot_ispec(s.dspec, s.ispec, s.uspec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=solve(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = source_init(s);
s = field_init(s);

s.it = 1;

s = kspec_init(s);
s.ispec = ispec_init(s.ispec, s); 
s.uspec = ispec_init(s.uspec, s); 
s.dspec = ispec_init(s.dspec, s); 

if 0, 
	plot_fields(s);
else 
	plot_kspec(s);
end

for it = 2:s.maxiter,

	s.it = it;

	[sn, sE] = source_update(s);

	n1 = s.n.*s.soundDecay1 + s.v.*s.soundDecay2 + sn;
	s.v = s.v.*s.soundDecay4 - s.h/2*s.k.^2.*(s.n.*s.soundDecay4 + n1);
	s.n = n1;

	s.E = s.E.*s.langmuirDecay + sE;

	s = kspec_update(s);
	s.ispec = ispec_update(s.ispec, s);
	s.uspec = ispec_update(s.uspec, s);
	s.dspec = ispec_update(s.dspec, s);

	if mod(s.it,100) == 0 | s.it == s.maxiter,
		fprintf(1,'Iter = %6d/%6d%s', s.it, s.maxiter, char(8*ones(1,20)));
		if 0, 
			plot_fields(s);
		else, 
			plot_kspec(s); 
		end
		drawnow
	end
end
fprintf(1,'\n');

s = kspec_normalise(s);
s.ispec = ispec_normalise(s.ispec, s);
s.uspec = ispec_normalise(s.uspec, s);
s.dspec = ispec_normalise(s.dspec, s);

plot_ispec(s.dspec, s.ispec, s.uspec);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = variable_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% physical constants
si=SIconstants;
s.me = si.me;
s.mi = s.miamu*si.amu;
s.Ne = s.ne*s.L_r^3;

% plasma dimensions
s.we = sqrt(s.ne*si.e^2/(si.eps0*s.me));
s.ve = sqrt(si.kb*s.Te/s.me);
s.eta = (s.Te+3*s.Ti)/s.Te;
s.lambdae = sqrt(si.eps0*si.kb*s.Te/(s.ne*si.e^2));

% normalisation constants
s.tau = 3/2*(s.mi/(s.me*s.eta))/s.we;                                       % [s]
s.chi = 3/2*sqrt(s.mi/(s.me*s.eta))*s.lambdae;                              % [m]
% TO CHECK !
s.epsilon = 4/sqrt(3)*s.eta*sqrt(s.me/s.mi)*sqrt(s.ne*si.kb*s.Te/si.eps0);  % [V m-1]
%s.epsilon = 2*sqrt(s.eta*s.me/s.mi)*sqrt(s.ne*si.kb*s.Te/si.eps0);  % [V m-1]
s.nu = 4/3*s.eta*(s.me/s.mi)*s.ne;                                          % [m-3]

s.k_r = 2*pi* [-s.N/2:s.N/2-1]/s.L_r;
s.M = floor(s.N/3);
s.pad0 = find([-s.N/2:s.N/2-1] > s.M | [-s.N/2:s.N/2-1] < -s.M);
s.h = s.h_r/s.tau;
s.k = s.k_r*s.chi;

s.k0 = find(s.k == 0);
s.ip = find(s.k > 0 & [-s.N/2:s.N/2-1] <= s.M);
s.im = find(s.k < 0 & [-s.N/2:s.N/2-1] >= -s.M);

s.ikc = find([-s.N/2:s.N/2-1]' >= -s.M & [-s.N/2:s.N/2-1]' <= s.M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = damping_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s.nui = zeros(size(s.k));
s.nui = s.nuic/2*s.tau + sqrt(pi/8)* ...
	(sqrt(s.me/s.mi)+(s.Te/s.Ti)^1.5*exp(-s.Te/(2*s.Ti)-3/2))*abs(s.k);

if 0
si=SIconstants;
w = s.k_r*sqrt(si.kb*(s.Te+3*s.Ti)/s.mi) - i*s.nui/s.tau;
for j=1:s.N,
	w(j) = disproot(s.k_r(j),w(j),s.ne,s.Te,s.Ti,s.miamu);
	s.nui(j) = s.nuic/2*s.tau - imag(w(i))*s.tau;
	fprintf(1,'.');
end
end

s.nue = zeros(size(s.k));
ik = find(s.k~=0);
s.nue(ik) = s.nuec/2*s.tau + sqrt(pi/8)*(3/2)^4*(s.mi/(s.eta*s.me))^2.5./ ...
	abs(s.k(ik)).^3.*(1+4/3*s.k(ik).^2*s.eta*s.me/s.mi).* ...
	exp(-9/8*s.mi/(s.eta*s.me)./s.k(ik).^2-3/2);

if 0
si=SIconstants;
w = s.we + 3/2*s.k_r*s.lambdae - i*s.nue/s.tau;
for j=1:s.N,
  w(j) = disproot(s.k_r(j),w(j),s.ne,s.Te,s.Ti,s.miamu);
	  s.nui(i) = s.nuic/2*s.tau - imag(w(j))*s.tau;
		  fprintf(1,'.');
end
end



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
a2 = (s.k_r*s.lambdae).^2;
%s.sn = sqrt(s.ne^2/s.Ne)*sqrt(1./(1+a2)./(1+a2+s.Te/s.Ti))/s.nu;
s.sn = sqrt((s.ne./s.nu)*(1./(1+a2)./(1+a2+s.Te/s.Ti)));
if 1
n1=sqrt(s.ne^2/s.Ne)*sqrt(1./(1+a2)./(1+a2+s.Te/s.Ti))/s.nu;
n2=sqrt((s.ne/s.nu)*(1./(1+a2)./(1+a2+s.Te/s.Ti)))*s.N;
n2([s.pad0 s.k0]) = 0.0;
u1=rand(size(s.im));
rp = zeros(size(s.k));
rp(s.im) = exp(i*2*pi*u1);
rp(s.ip) = flipud(conj(rp(s.im)));
subplot(211), plot(s.k_r*s.lambdae, n2.^2/s.ne); 
subplot(212), plot(real(ifft(n2.*rp)));
pause
end
% multiply by s.N for compatibility with the definition of inverse DFT
s.sn = s.sn*s.N;
s.sn([s.pad0 s.k0]) = 0.0;
s.Dn = sqrt(2*s.nui.*s.k.^2./(s.k.^2+4*s.nui.^2)).*s.sn;
s.Dn([s.pad0 s.k0]) = 0.0;
%subplot(111), plot(s.k, [s.sn; s.Dn]); pause

%s.sE = sqrt(s.ne^2/s.Ne)*(si.e/si.eps0)*s.lambdae./sqrt(2*(1+a2));
s.sE = sqrt((s.ne/s.nu)*(si.e/si.eps0)^2*(s.lambdae/s.chi)^2./(2*(1+a2)));
if 1
E1=sqrt(s.ne^2/s.Ne)*(si.e/si.eps0)*(s.lambdae/s.chi)./sqrt(2*(1+a2));
E2=sqrt((s.ne/s.nu)*(si.e/si.eps0)^2*(s.lambdae/s.chi)^2./(2*(1+a2)));
E2([s.pad0 s.k0]) = 0.0;
subplot(212), plot(s.k_r*s.lambdae, (abs(s.k).*E2*si.eps0/si.e).^2/s.ne)
pause
end
s.sE = s.sE/s.epsilon*s.N;
s.sE([s.pad0 s.k0]) = 0.0;
s.DE = sqrt(s.nue).*s.sE;
s.DE([s.pad0 s.k0]) = 0.0;
%subplot(111), plot(s.k, [s.sE; s.DE]); pause

rand('twister', s.seeds(s.is));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sn, sE] = source_update(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u1 = rand(size(s.im));
u2 = rand(size(s.im));
Wh = zeros(size(s.k));
if 0
Wh(s.im) = sqrt(s.h)*sqrt(-2*log(u1)).*exp(i*2*pi*u2);
Wh(s.ip) = flipud(conj(Wh(s.im)));
sn = s.Dn.*Wh;
else
Wh(s.im) = sqrt(-2*log(u1)).*exp(i*2*pi*u2);
Wh(s.ip) = flipud(conj(Wh(s.im)));
sn = sqrt(s.Dn.^2./(2*s.nui).*(1-exp(-2*s.nui*s.h))).*Wh;
end

u1 = rand(size(s.k));
u2 = rand(size(s.k));
if 0
Wh = sqrt(s.h)*sqrt(-2*log(u1)).*exp(i*2*pi*u2);
sE = s.DE.*Wh;
else
Wh = sqrt(-2*log(u1)).*exp(i*2*pi*u2);
sE = sqrt(s.DE.^2./(2*s.nue).*(1-exp(-2*s.nue*s.h))).*Wh;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = field_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%

s.theta = acos(-sqrt(4*s.nui.^2./(s.k.^2+4*s.nui.^2)));
s.n = s.sn.*exp(i*s.theta/2);
s.v = s.Dn.*abs(s.k)./sqrt(2*s.nui).*exp(-i*s.theta/2);
s.E = s.sE;

%plot(s.k, 180/pi*acos(4*sign(s.k).*s.nui.^2./(s.k.^2+4*s.nui.^2)))
%pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = kspec_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s.spec.k = s.k_r(s.ikc);
if 0
s.spec.nk2 = real(s.n(s.ikc).*conj(s.n(s.ikc)));
s.spec.Ek2 = real(s.E(s.ikc).*conj(s.E(s.ikc)));
s.spec.nspec = 1;
else
s.spec.nk2 = zeros(length(s.spec.k),s.maxiter,s.ns);
s.spec.Ek2 = zeros(length(s.spec.k),s.maxiter,s.ns);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = kspec_update(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
s.spec.nk2 = s.spec.nk2 + real(s.n.*conj(s.n));
s.spec.Ek2 = s.spec.Ek2 + real(s.E.*conj(s.E));
s.spec.nspec = s.spec.nspec+1;
else
s.spec.nk2(:,s.it,s.is) = real(s.n(s.ikc).*conj(s.n(s.ikc)));
s.spec.Ek2(:,s.it,s.is) = real(s.E(s.ikc).*conj(s.E(s.ikc)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = kspec_normalise(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
s.spec.nk2 = s.spec.nk2/s.N^2*s.nu^2;
else
s.spec.nk2 = s.spec.nk2/s.N*s.nu;
end
s.spec.Ek2 = s.spec.Ek2/s.N^2*s.epsilon^2;

s.sn = s.sn/s.N*s.nu;
s.sE = s.sE/s.N*s.epsilon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spec = ispec_init(spec,s) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec.k = s.k_r(spec.ik);
spec.ts = spec.stride*s.h_r;
spec.fs = 1/spec.ts;
spec.f = spec.fs*[-spec.ns/2:spec.ns/2-1]/spec.ns;
spec.skf2 = zeros(length(spec.k), length(spec.f));
spec.skt = zeros(length(spec.k), length(spec.f));
spec.nspec = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spec = ispec_update(spec,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rem(s.it, spec.stride) == 0,
	it = rem(fix(s.it/spec.stride),spec.ns)+1;
	if strcmp(spec.type,'i'),
		spec.skt(:, it) = s.n(spec.ik);
	else,
		spec.skt(:, it) = s.E(spec.ik);
	end
	if it == spec.ns,
		if strcmp(spec.type,'d'),
			skf = fft(spec.skt.').';
			spec.skf2 = spec.skf2 + real(skf.*conj(skf));
		else
			skf = spec.ns * ifft(spec.skt.').';
			spec.skf2 = spec.skf2 + real(skf.*conj(skf));
		end
		spec.nspec = spec.nspec+1;
		[f,k] = meshgrid(spec.f, spec.k);
		spec.sk2 = integrate(f', spec.skf2');
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spec = ispec_normalise(spec,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if spec.nspec ~= 0,
	spec.skf2 = spec.skf2/s.N^2/spec.ns^2/spec.nspec;
end
if strcmp(spec.type,'i'),
	spec.skf2 = fftshift(spec.skf2*s.nu^2, 2);
else
	k2 = repmat(s.k_r(spec.ik).^2',1,spec.ns);
	spec.skf2 = fftshift(k2.*spec.skf2*s.epsilon^2, 2);
	si = SIconstants;
	spec.skf2 = spec.skf2*s.me/s.mi*(si.eps0/(si.kb*s.Te))^2;
end
[f,k] = meshgrid(spec.f, spec.k);
spec.sk2 = integrate(f', spec.skf2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_ispec(dspec,ispec,uspec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)

subplot(231);
imagesc(dspec.f, dspec.k, dspec.skf2); 
colorbar

subplot(232); 
imagesc(ispec.f, ispec.k, ispec.skf2);
axis xy
colorbar

subplot(233);
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
figure(1)

subplot(211), 
plot(s.k, real(s.n.*conj(s.n))*(s.nu/s.N)^2);
subplot(212), 
plot(s.k, real(s.E.*conj(s.E))*(s.epsilon/s.N)^2);

%%%%%%%%%%%%%%%%%%%%%%
function plot_kspec(s)
%%%%%%%%%%%%%%%%%%%%%%
figure(1)

if size(s.spec.nk2,2)==0 || size(s.spec.Ek2,2)==0, 
  return;
end

subplot(211),
h=plot(s.spec.k, mean(s.spec.nk2(:,1:s.it,s.is),2)*(s.nu/s.N)^2, ...
       s.spec.k,s.sn(s.ikc).^2*(s.nu/s.N)^2);
set(h(2),'linewidth',2)
set(gca,'xlim',[min(s.spec.k) max(s.spec.k)]);
title(sprintf('%d paths @ h=%.2g T=%.2g', s.ns, s.h, s.h*s.it))
ylabel('<|n_k(t)|^2>_t')
set(gca,'xticklabel',[]);

subplot(212),
h=plot(s.spec.k, mean(s.spec.Ek2(:,1:s.it,s.is),2)*(s.epsilon/s.N)^2, ...
       s.spec.k,s.sE(s.ikc).^2*(s.epsilon/s.N)^2);
set(h(2),'linewidth',2)
set(gca,'xlim',[min(s.spec.k) max(s.spec.k)]);
ylabel('<|E_k(t)|^2>_t')
xlabel('k')

