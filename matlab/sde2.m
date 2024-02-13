function s = sde2(s)
% function s = sde2(s)

%
% $Id: sde2.m,v 1.11 2011/03/26 09:20:42 patrick Exp $
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

if nargin == 0 | (nargin == 1 & ~isfield(s,'h')),
	s.h = 1e-2/s.c;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'maxiter')),
	s.maxiter = 2500*s.c;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'gamma')),
	s.gamma = 1.0;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'beta')),
	s.beta = sqrt(4*pi*s.gamma/10);
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'phi0')),
	s.phi0 = 0;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'dphi0')),
	s.dphi0 = 0;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'rng')),
	s.rng='normal';
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'seeds')),
	s.seeds = 1:1000;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'k')),
	s.k = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1 & isfield(s,'bench'),
	s.seed = 1;
	s = bench_sde2(s);
	return
end

s.phit2m = zeros(s.maxiter,1);
s.dphit2m = zeros(s.maxiter,1);
s.nm = 0;

fprintf(1,'seeds =');
for seed = s.seeds,
	s.seed = seed;
	s = solve_sde2(s);
	s.phit2m = s.phit2m + s.phit2;
	s.dphit2m = s.dphit2m + s.dphit2;
	s.nm = s.nm + 1;
	if rem(seed,100)==0,
		fprintf(1,' %d', seed);
	end
end
fprintf(1,'\n');

s.phit2m = s.phit2m/s.nm;
s.dphit2m = s.dphit2m/s.nm;

s = analytic_init(s);

if 0
	subplot(311), 
	plot(s.t,s.phit2)
	ylabel('|\phi(t)|^2');
	title(sprintf('sde2: %d paths @ %d points (c=%.1g), rng=%s, k=%.2g', ...
		length(s.seeds),length(s.phit2m), s.c, s.rng, s.k))
	subplot(312), 
	plot(s.t,s.phit2m)
	ylabel('<|\phi(t)|^2>');
	subplot(313),
	if 1
		plot(s.t,s.dphit2m,s.t,s.dphit2a)
		ylabel('<|d\phi(t)|^2>');
	else
		plot(s.t,s.phit2a)
		ylabel('<|\phi(t)|^2>');
	end
	xlabel('t')
else
	subplot(211), 
	plot(s.t,s.phit2m,s.ta,s.phit2a,'linewidth',1)
	title(sprintf('sde2: %d paths @ %d points (c=%.1g), rng=%s, k=%.2g', ...
		length(s.seeds),length(s.phit2m), s.c, s.rng, s.k));
	ylabel('<|\phi(t)|^2>');
	subplot(212), 
	plot(s.t,s.dphit2m,s.ta,s.dphit2a,'linewidth',1)
	ylabel('<d|\phi(t)|^2>');
	xlabel('t')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = solve_sde2(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%
frq = s.k.^2-s.gamma.^2;
s.ct = ones(size(s.k));
s.st = repmat(s.h, size(s.k));
ip = find(frq>0);
s.ct(ip) = cos(sqrt(frq(ip))*s.h);
s.st(ip) = sin(sqrt(frq(ip))*s.h)./sqrt(frq(ip));
im=find(frq<0);
s.ct(im) = cosh(sqrt(-frq(im))*s.h);
s.st(im) = sinh(sqrt(-frq(im))*s.h)./sqrt(-frq(im));

s.nu1 = exp(-s.gamma*s.h).*(s.ct+s.gamma.*s.st);
s.nu2 = exp(-s.gamma*s.h).*s.st;
s.nu4 = exp(-2*s.gamma*s.h);

s = source_init(s);
s.t = [0:s.maxiter-1]'*s.h;
s.phi = zeros(s.maxiter,1);
s.dphi = zeros(s.maxiter,1);

if 0,
	s.phi(1) = s.phi0;
	s.dphi(1) = s.dphi0;
	s = method1(s);
%	phi1 = s.phi.*conj(s.phi);
%	dphi1 = s.dphi.*conj(s.dphi);
end

if 1,
	if ~isfield(s,'A'),
		s.A = spdiags([[repmat([-1;-1+2*s.h*s.gamma],s.maxiter-1,1);0;0],...
								 	[0;repmat(s.h*[-1;s.k^2],s.maxiter-1,1);0],...
								 	[ones(2*s.maxiter,1)]],...
								 	[-2:0], 2*s.maxiter, 2*s.maxiter);
		[s.L,s.U]=luinc(s.A,1e-6);
	end
	x = zeros(2*(s.maxiter-1),1);
	x(1:2:end-1) = s.source1(1:end-1);
	s.b = [s.phi0; s.dphi0; x];
	s = method3(s);
%	phi2 = s.phi.*conj(s.phi);
%	dphi2 = s.dphi.*conj(s.dphi);
end
%subplot(211), plot(s.t,[phi1,phi2]);
%subplot(212), plot(s.t,[dphi1,dphi2]);
%pause

s.phit2 = s.phi.*conj(s.phi);
s.dphit2 = s.dphi.*conj(s.dphi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = source_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rand('state', s.seed);
for n=1:2,
	switch s.rng
		case 'rpa',
			x1 = rand(s.maxiter,1);
			x2 = rand(s.maxiter,1);
			xi = exp(i*2*pi*x2);
		case 'normal',
			x1 = rand(s.maxiter,1);
			x2 = rand(s.maxiter,1);
			xi = sqrt(-2*log(x1)).*exp(i*2*pi*x2);
			xi = xi/sqrt(2);
	end
	if n==1,
		s.source1 = sqrt(s.h).*s.beta.*xi;
	elseif n==2,
		s.source2 = sqrt(s.h).*s.beta.*xi;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = analytic_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global nui k sign2 sigv2
nui=s.gamma;
sign2=s.beta^2;
sigv2=0;
k=s.k;
y0=[s.phi0*conj(s.phi0);s.dphi0*conj(s.dphi0); ...
		s.phi0*conj(s.dphi0)+conj(s.phi0)*s.dphi0];
if 0
	[t,y]=ode45('ode1',[min(s.t) max(s.t)],y0);
else
	t = s.t;
	y=ode1a(t,y0);
end
s.ta = t;
s.phit2a = y(:,1);
s.dphit2a = y(:,2);


function s = method1(s)
for j = 1:s.maxiter-1, 
	s.phi(j+1) = s.phi(j).*s.nu1+s.dphi(j)*s.nu2+s.source1(j);
	s.dphi(j+1) = s.dphi(j)*s.nu4-(s.h/2)*s.k.^2.*(s.phi(j).*s.nu4+s.phi(j+1));
end
phi1 = s.phi.*conj(s.phi);
dphi1 = s.dphi.*conj(s.dphi);
for j = 1:s.maxiter-1,
	s.phi(j+1) = s.phi(j)+s.dphi(j)*s.h+s.source1(j);
	s.dphi(j+1) = (1-2*s.gamma*s.h)*s.dphi(j)-s.k^2*s.phi(j)*s.h;
end
phi2 = s.phi.*conj(s.phi);
dphi2 = s.dphi.*conj(s.dphi);
subplot(211), plot(s.t,[phi1 phi2]);
subplot(212), plot(s.t,[dphi1 dphi2]);
pause

function s = method2(s)
R = qr(s.A);
s.u = R\(R'\(s.A'*s.b));
r = s.b-s.A*s.u;
e = R\(R'\(s.A'*r));
s.u = s.u+e;
s.phi = s.u(1:2:end-1);
s.dphi = s.u(2:2:end);

function s = method3(s)
[s.u,flag,relres,iter] = bicgstab(s.A,s.b,1e-6,10,s.L,s.U);
conv_check(flag,relres,iter);
s.phi = s.u(1:2:end-1);
s.dphi = s.u(2:2:end);

function s = method4(s)
[s.u,flag,relres,iter] = bicg(s.A,s.b,1e-4,10,s.L,s.U);
conv_check(flag,relres,iter);
s.phi = s.u(1:2:end-1);
s.dphi = s.u(2:2:end);

function s = method5(s)
[s.u,flag,relres,iter] = cgs(s.A,s.b,1e-6,10,s.L,s.U);
conv_check(flag,relres,iter),
s.phi = s.u(1:2:end-1);
s.dphi = s.u(2:2:end);

function s = method6(s)
[s.u,flag,relres,iter] = gmres(s.A,s.b,5,1e-6,10,s.L,s.U);
conv_check(flag,relres,iter);
s.phi = s.u(1:2:end-1);
s.dphi = s.u(2:2:end);

function s = method7(s)
[s.u,flag,relres,iter] = qmr(s.A,s.b,1e-6,10,s.L,s.U);
conv_check(flag,relres,iter);
s.phi = s.u(1:2:end-1);
s.dphi = s.u(2:2:end);

function conv_check(flag,relres,iter)
if flag ~= 0,
	error(sprintf('flag=%d,relres=%e,iter=%d',flag,relres,iter));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = bench_sde2(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%
s.nu = exp(-s.gamma*s.h);
s = source_init(s);
s.t = [0:s.maxiter-1]'*s.h;
s.phi = zeros(s.maxiter,1);
s.dphi = zeros(s.maxiter,1);

s.nloop = 1;
s.nmethod = 7;
s.PHIT2 = zeros(s.maxiter,s.nmethod);

% method 1
s.phi(1) = s.phi0;
s.dphi(1) = s.dphi0;

% method 2-7
s.A = spdiags([[repmat([-1;-1+2*s.h*s.gamma],s.maxiter-1,1);0;0],...
							 [0;repmat(s.h*[-1;s.k^2],s.maxiter-1,1);0],...
							 [ones(2*s.maxiter,1)]],...
							 [-2:0], 2*s.maxiter, 2*s.maxiter);
[s.L,s.U]=luinc(s.A,1e-6);
x = zeros(2*(s.maxiter-1),1);
x(1:2:end-1) = s.source1(1:end-1);
s.b = [s.phi0; s.dphi0; x];

for j=1:s.nmethod,
	tic;
	eval(['for i=1:s.nloop, s = method' num2str(j) '(s); end;']);
	t(j) = toc;
	fprintf(1,'method %d: %.2f s\n', j, t(j));
	s.PHIT2(:,j) = s.phi.*conj(s.phi);
end

subplot(211)
plot(s.t,s.PHIT2)
xlabel('t');
ylabel('|\phi(t)|^2')
title(sprintf('sde2: %d loops @ %d points, k=%.2g', s.nloop, s.maxiter, s.k))
for i=1:s.nmethod,
	leg{i} = sprintf('method %d (%.0f ms)', i, 100*t(i)/s.nloop);
end
legend(leg,1)
orient tall

subplot(212)
semilogy(s.t,abs(s.PHIT2(:,:)-repmat(s.PHIT2(:,1),1,s.nmethod)))
ylabel('||\phi_i(t)|^2-|\phi_1(t)|^2|')
xlabel('t');
set(gca,'ylim',[1e-30 1e-3]);
clear leg
for i=1:6,
	leg{i} = sprintf('method %d - method 1', i+1);
end
legend(leg,4)
orient tall
print -dpsc sde2-bench

function s=sparse2_init(s)
s.A = spdiags([[repmat(s.h*s.k^2,s.maxiter-1,1);zeros(s.maxiter+1,1)],...
               [repmat(-s.h,s.maxiter-1,1);0;...
							 repmat(-1+2*s.h*s.gamma,s.maxiter-1,1);0],...
               [ones(2*s.maxiter,1)],...
							 [zeros(s.maxiter+1,1);repmat(-s.h,s.maxiter-1,1)]],...
               [-s.maxiter-1 -1:0 s.maxiter], 2*s.maxiter, 2*s.maxiter);
s.b = [s.phi0; s.source1(1:end-1); s.dphi0; zeros(s.maxiter-1,1)];

