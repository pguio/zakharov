function s = sde0(s)
% function s = sde0(s)

%
% $Id: sde0.m,v 1.16 2011/03/26 09:20:41 patrick Exp $
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

if nargin == 0 | (nargin == 1 & ~isfield(s,'k')),
	s.k = 1.0;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'beta')),
	s.beta = sqrt(4*pi*s.gamma/10);
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'phi0')),
	s.phi0 = 0;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'rng')),
	s.rng='normal';
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'seeds')),
	s.seeds = 1:1000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1 & isfield(s,'bench'),
	s.seed = 1;
	s = bench_sde0(s);
	return
end

s.phit2m = zeros(s.maxiter,1);
s.nphit2 = 0;

fprintf(1,'seeds =');
for seed = s.seeds,
	s.seed = seed;
	s = solve_sde0(s);
	s.phit2m = s.phit2m + s.phit2;
	s.nphit2 = s.nphit2 + 1;
	if rem(seed,100)==0,
		fprintf(1,' %d', seed);
	end
end
fprintf(1,'\n');

s.phit2m = s.phit2m/s.nphit2;

s = analytic_init(s);

if 0
	subplot(311), 
	plot(s.t,s.phit2)
	ylabel('|\phi(t)|^2');
	title(sprintf('sde0: %d paths @ %d points (c=%.1g), rng=%s', ...
		length(s.seeds),length(s.phit2m), s.c, s.rng))
	subplot(312), 
	plot(s.t,s.phit2m)
	ylabel('<|\phi(t)|^2>');
	subplot(313),
	plot(s.t,s.phit2a)
	xlabel('t')
	ylabel('<|\phi(t)|^2>');
else
	subplot(111),
	plot(s.t,s.phit2m,s.t,s.phit2a,'linewidth',1)
	xlabel('t')
	ylabel('|<\phi(t)>|^2');
	title(sprintf('sde0: %d paths @ %d points (c=%.1g), rng=%s', ...
		length(s.seeds),length(s.phit2m), s.c, s.rng))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = solve_sde0(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%
s.nu = exp(-(i*s.k^2+s.gamma)*s.h);
s = source_init(s);
s.t = [0:s.maxiter-1]'*s.h;
s.phi = zeros(s.maxiter,1);

if 0, % yields method 1
s.phi(1) = s.phi0;
s = method1(s);
end

if 0, % yields method 2 to 7
if ~isfield(s,'A'),
	s.A = spdiags([[repmat(-s.nu,s.maxiter-1,1);0],[ones(s.maxiter,1)]],...
                [-1:0], s.maxiter, s.maxiter);
	[s.L,s.U]=luinc(s.A,1e-6);
end
s.b = [s.phi0; s.source(1:end-1)];
s = method3(s);
end

if 1, % yields method 8
s.b = [s.phi0; s.source(1:end-1)];
s = method8(s);
end

s.phit2 = s.phi.*conj(s.phi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = source_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
rand('state', s.seed);
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
if 0
	s.source = sqrt(s.h).*s.beta.*xi;
else
	s.source = sqrt(s.h).*s.beta.*xi* ...
		(1+(1+i*s.k^2+s.gamma)*exp(-(i*s.k^2+s.gamma)*s.h));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = analytic_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s.phit2a = exp(-2*s.gamma*s.t).*s.phi0.^2 + ...
          s.beta.^2./2/s.gamma.*(1-exp(-2*s.gamma*s.t));


function s = method1(s)
for j = 1:s.maxiter-1, 
	s.phi(j+1) = s.phi(j).*s.nu + s.source(j);
end
subplot(211), plot(s.t, s.phi.*conj(s.phi));
for j = 1:s.maxiter-1, 
	s.phi(j+1) = s.phi(j).*(1-(i*s.k^2+s.gamma)*s.h) + s.source(j);
end
subplot(212), plot(s.t, s.phi.*conj(s.phi)); pause

function s = method2(s)
R = qr(s.A);
s.phi = R\(R'\(s.A'*s.b));
r = s.b-s.A*s.phi;
e = R\(R'\(s.A'*r));
s.phi = s.phi+e;

function s = method3(s)
[s.phi,flag,relres,iter] = bicgstab(s.A,s.b,1e-6,10,s.L,s.U);
conv_check(flag,relres,iter);

function s = method4(s)
[s.phi,flag,relres,iter] = bicg(s.A,s.b,1e-6,10,s.L,s.U);
conv_check(flag,relres,iter);

function s = method5(s)
[s.phi,flag,relres,iter] = cgs(s.A,s.b,1e-6,10,s.L,s.U);
conv_check(flag,relres,iter),

function s = method6(s)
[s.phi,flag,relres,iter] = gmres(s.A,s.b,5,1e-6,10,s.L,s.U);
conv_check(flag,relres,iter);

function s = method7(s)
[s.phi,flag,relres,iter] = qmr(s.A,s.b,1e-6,10,s.L,s.U);
conv_check(flag,relres,iter);

function s = method8(s)
s.phi = filter(1,[1 -s.nu],s.b);

function conv_check(flag,relres,iter)
if flag ~= 0,
	error(sprintf('flag=%d,relres=%e,iter=%d',flag,relres,iter));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = bench_sde0(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%
s.nu = exp(-s.gamma*s.h);
s = source_init(s);
s.t = [0:s.maxiter-1]'*s.h;
s.phi = zeros(s.maxiter,1);

s.nloop = 100;
s.nmethod = 8;
s.PHIT2 = zeros(s.maxiter,s.nmethod);

% method 1
s.phi(1) = s.phi0;

% method 2-7
s.A = spdiags([[repmat(-s.nu,s.maxiter-1,1);0],[repmat(1,s.maxiter,1)]],...
              [-1:0], s.maxiter, s.maxiter);
[s.L,s.U] = luinc(s.A,1e-6);

% method 2-8
s.b = [s.phi0; s.source(1:end-1)];

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
title(sprintf('sde0: %d loops @ %d points', s.nloop, s.maxiter))
for i=1:s.nmethod,
	leg{i} = sprintf('method %d (%.0f ms)', i, 1000*t(i)/s.nloop);
end
legend(leg,1)
orient tall

subplot(212)
semilogy(s.t,abs(s.PHIT2(:,:)-repmat(s.PHIT2(:,1),1,s.nmethod)))
ylabel('||\phi_i(t)|^2-|\phi_1(t)|^2|')
xlabel('t');
set(gca,'ylim',[1e-30 1e-3]);
clear leg
for i=1:s.nmethod-1,
	leg{i} = sprintf('method %d - method 1', i+1);
end
legend(leg,4)
orient tall
print -dpsc sde0-bench

% plot(s.t,abs(s.PHIT2(:,8)-s.PHIT2(:,1)))
