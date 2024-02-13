function s = sde6(s)
% function s = sde6(s)

%
% $Id: sde6.m,v 1.20 2023/02/23 08:57:34 patrick Exp $
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

if nargin == 0 | (nargin == 1 & ~isfield(s,'h')),
	s.h = 1e-2;
end
if nargin == 0 | (nargin == 1 & ~isfield(s,'maxiter')),
	s.maxiter = 2500;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'seeds')),
	s.seeds = 1:300;
end

if nargin == 0 | (nargin == 1 & ~isfield(s,'nk')),
	s.nk = 4;
end
if nargin == 0 | (nargin == 1 & ~isfield(s,'nk')),
  s.k = [1 10 20 50];
end
if nargin == 0 | (nargin == 1 & ~isfield(s,'nu')),
	s.nu = [1 10 100 200];
end
if nargin == 0 | (nargin == 1 & ~isfield(s,'phi0')),
	s.phi0 = [1 1 1 1];
end
if nargin == 0 | (nargin == 1 & ~isfield(s,'dphi0')),
	s.dphi0 = [1 1 1 1];
end
if nargin == 0 | (nargin == 1 & ~isfield(s,'D')),
  s.D = { sqrt([1 10 100 200]), [1 1 1 1] };
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,' SDE for ion-acoustic waves\n');
fprintf(1,' dn =       v dt            + Dn dW\n');
fprintf(1,' dv = -2 nu v dt - k^2 n dt + Dv dW\n\n');
fprintf(1,' nu    = '); fprintf(1,'%10.4g ', s.nu); fprintf(1,'\n');
fprintf(1,' k     = '); fprintf(1,'%10.4g ', s.k); fprintf(1,'\n');
fprintf(1,' Dn    = '); fprintf(1,'%10.4g ', s.D{1}); fprintf(1,'\n');
fprintf(1,' Dv    = '); fprintf(1,'%10.4g ', s.D{2}); fprintf(1,'\n');
fprintf(1,' Re n0 = '); fprintf(1,'%10.4g ', real(s.phi0)); fprintf(1,'\n');
fprintf(1,' Im n0 = '); fprintf(1,'%10.4g ', imag(s.phi0)); fprintf(1,'\n');
fprintf(1,' Re v0 = '); fprintf(1,'%10.4g ', real(s.dphi0)); fprintf(1,'\n');
fprintf(1,' Im v0 = '); fprintf(1,'%10.4g ', imag(s.dphi0)); fprintf(1,'\n\n');


s.ns = length(s.seeds);
s.phi2s = zeros(s.maxiter, s.nk, s.ns);
s.dphi2s = zeros(s.maxiter, s.nk, s.ns);
s.phidphi2s = zeros(s.maxiter, s.nk, s.ns);
s.t = [0:s.maxiter-1]'*s.h;
moduloDisplay = fix(s.ns/10);
fprintf(1,' seed = ');
for seed = s.seeds,
	s.seed = seed;
	s = solve_sde6(s);
	s.phi2s(:,:,seed) = s.phi2;
	s.dphi2s(:,:,seed) = s.dphi2;
	s.phidphi2s(:,:,seed) = s.phi.*conj(s.dphi)+conj(s.phi).*s.dphi;
	if ~rem(seed,moduloDisplay),
		fprintf(1,'%4d/%4d%s', seed, s.ns, char(8*ones(1,9)));
	end
end
fprintf(1,'\n\n');

s = analytic_init(s);

s.phi2stats = zeros(s.maxiter,s.nk, 6);
s.phi2stats(:,:,1) = mean(s.phi2s,3);
stdstats = std(s.phi2s,0,3)/sqrt(2);
s.phi2stats(:,:,2) = s.phi2stats(:,:,1) + stdstats;
s.phi2stats(:,:,3) = s.phi2stats(:,:,1) - stdstats;
s.phi2stats(:,:,4:6) = quantile(s.phi2s,[0.25 0.50 0.75],3);

s.dphi2stats = zeros(s.maxiter,s.nk, 6);
s.dphi2stats(:,:,1) = mean(s.dphi2s,3);
stdstats = std(s.dphi2s,0,3)/sqrt(2);
s.dphi2stats(:,:,2) = s.dphi2stats(:,:,1) + stdstats;
s.dphi2stats(:,:,3) = s.dphi2stats(:,:,1) - stdstats;
s.dphi2stats(:,:,4:6) = quantile(s.dphi2s,[0.25 0.50 0.75],3);

s.phidphi2stats = zeros(s.maxiter,s.nk, 6);
s.phidphi2stats(:,:,1) = mean(s.phidphi2s,3);
stdstats = std(s.phidphi2s,0,3)/sqrt(2);
s.phidphi2stats(:,:,2) = s.phidphi2stats(:,:,1) + stdstats;
s.phidphi2stats(:,:,3) = s.phidphi2stats(:,:,1) - stdstats;
s.phidphi2stats(:,:,4:6) = quantile(s.phidphi2s,[0.25 0.50 0.75],3);


s=plotsde(s);



%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = solve_sde6(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise \Delta W(t,k) 
rand('twister', s.seed);
for n=1:2,
  u1 = rand(s.maxiter,s.nk);
	u2 = rand(s.maxiter,s.nk);
	if 0
	% Euler-Maruyama first order scheme
	s.wh{n} = sqrt(s.h)*sqrt(-2*log(u1)).*exp(i*2*pi*u2);
	s.source{n} = repmat(s.D{n}, s.maxiter,1).*s.wh{n};
	else
	% exact scheme from Gillespie (1997)
	s.wh{n} = sqrt(-2*log(u1)).*exp(i*2*pi*u2);
	s.source{n} = repmat(sqrt(s.D{n}.^2./s.nu/2.*(1-exp(-2*s.nu*s.h))), s.maxiter,1).*s.wh{n};
	end
end

s = solver_init(s);

% allocate phi(t,k)
s.phi = zeros(s.maxiter, s.nk);
s.dphi = zeros(s.maxiter, s.nk);

if 0
% explicit method
	s.phi(1,:) = s.phi0;
	s.dphi0(1,:) = s.dphi0;
	for j = 2:s.maxiter,
		s.phi(j,:) = s.phi(j-1,:).*s.nu1+s.dphi(j-1,:).*s.nu2+s.source{1}(j,:);
		s.dphi(j,:) = s.dphi(j-1,:).*s.nu4-...
			s.h/2*s.k.^2.*(s.phi(j-1,:).*s.nu4+s.phi(j,:))+s.source{2}(j,:);
	end
end

if 1
% slightly faster implicit method
	for k=1:s.nk,
	A = spdiags([[repmat([-1;-1+2*s.h*s.nu(k)],s.maxiter-1,1);0;0],...
              [0;repmat(s.h*[-1;s.k(k)^2],s.maxiter-1,1);0],...
	 					 [ones(2*s.maxiter,1)]],...
						 [-2:0], 2*s.maxiter, 2*s.maxiter);
	%[D,U]=luinc(A,1e-6);
	[D,U]=lu(A);
	x = zeros(2*(s.maxiter-1),1);
	x(1:2:end-1) = s.source{1}(2:end,k);
	x(2:2:end) = s.source{2}(2:end,k);
	b = [s.phi0(k); s.dphi0(k); x];
	[u,flag,relres,iter] = bicgstab(A,b,1e-6,10,D,U);
	s.phi(:,k) = u(1:2:end-1);
	s.dphi(:,k) = u(2:2:end);
	end
end

s.phi2 = s.phi.*conj(s.phi);
s.dphi2 = s.dphi.*conj(s.dphi);

%%%%%%%%%%%%&&%%%%%%%%%%%%%
function s = solver_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
frq = s.k.^2-s.nu.^2;
s.ct = ones(size(s.k));
s.st = repmat(s.h, size(s.k));
% be careful cos/cosh and sin/sinh!
ip = find(frq>0);
s.ct(ip) = cos(sqrt(frq(ip))*s.h);
s.st(ip) = sin(sqrt(frq(ip))*s.h)./sqrt(frq(ip));
im=find(frq<0);
s.ct(im) = cosh(sqrt(-frq(im))*s.h);
s.st(im) = sinh(sqrt(-frq(im))*s.h)./sqrt(-frq(im));

s.nu1 = exp(-s.nu*s.h).*(s.ct+s.nu.*s.st);
s.nu2 = exp(-s.nu*s.h).*s.st;
s.nu4 = exp(-2*s.nu*s.h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = analytic_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global nui k sign2 sigv2
for i=1:s.nk,
	k=s.k(i);
	y0=[s.phi0(i)*conj(s.phi0(i));s.dphi0(i)*conj(s.dphi0(i));  ...
		s.phi0(i)*conj(s.dphi0(i))+conj(s.phi0(i))*s.dphi0(i)];
	nui=s.nu(i);
	sign2=2*s.D{1}(i)^2;
	sigv2=2*s.D{2}(i)^2;
	if 0,
		[t,y]=ode45('ode1',[min(s.t) max(s.t)],y0);
	else
		t = s.t;
		y=ode1a(t,y0);
	end
	s.ta{i} = t;
	s.phi2a{i} = y(:,1);
	s.dphi2a{i} = y(:,2);
	s.phidphi2a{i} = y(:,3);
end

%%%%%%%%%%%%%%%%%%%%%
function s=plotsde(s)
%%%%%%%%%%%%%%%%%%%%%

if ~isfield(s,'nsub'),
  s.nsub = s.nk;
end

s.fig1=figure;
s.fig2=figure;
s.fig3=figure;
c=[0 0 1]; g=[0 .5 0]; r=[1 0 0];
for i=1:s.nk,
	%subplot(s.nk,3,3*i-2),
	figure(s.fig1)
	subplot(s.nsub,1,i)
	h=plot(s.t,squeeze(s.phidphi2stats(:,i,:)),s.ta{i},s.phidphi2a{i});
  set(h(1:3),'color',c);
  set(h(4:6),'color',g,'linestyle','--');
  set(h(7),'color',r,'linewidth',.6);
	set(gca,'xlim',[min(s.t) max(s.t)]);
	set(gca,'YAxisLocation','right');
	title(sprintf('\\nu=%.2g k=%.2g |n0|=%.2g |v_0|=%.2g D_n=%.2g D_v=%.2g ', ...
	  s.nu(i), s.k(i), abs(s.phi0(i)), abs(s.dphi0(i)), s.D{1}(i), s.D{2}(i)));
	if i==s.nk,
		xlabel('t');
	else
		set(gca,'xticklabel',[]);
	end
	%subplot(s.nk,3,3*i-1),
	figure(s.fig2)
	subplot(s.nsub,1,i)
	h=plot(s.t,squeeze(s.dphi2stats(:,i,:)),s.ta{i},s.dphi2a{i});
	set(h(1:3),'color',c);
	set(h(4:6),'color',g,'linestyle','--');
	set(h(7),'color',r,'linewidth',.6);
	set(gca,'xlim',[min(s.t) max(s.t)]);
	set(gca,'YAxisLocation','right');
	title(sprintf('\\nu=%.2g k=%.2g |n0|=%.2g |v_0|=%.2g D_n=%.2g D_v=%.2g ', ...
	    s.nu(i), s.k(i), abs(s.phi0(i)), abs(s.dphi0(i)), s.D{1}(i), s.D{2}(i)));
	if i==s.k,
		xlabel('t');
	else
		set(gca,'xticklabel',[]);
	end
	%subplot(s.nk,3,3*i),
	figure(s.fig3)
	subplot(s.nsub,1,i)
	h=plot(s.t,squeeze(s.phi2stats(:,i,:)),s.ta{i},s.phi2a{i});
	set(h(1:3),'color',c);
	set(h(4:6),'color',g,'linestyle','--');
	set(h(7),'color',r,'linewidth',.6);
	set(gca,'xlim',[min(s.t) max(s.t)]);
	set(gca,'YAxisLocation','right');
	title(sprintf('\\nu=%.2g k=%.2g |n0|=%.2g |v_0|=%.2g D_n=%.2g D_v=%.2g ', ...
	    s.nu(i), s.k(i), abs(s.phi0(i)), abs(s.dphi0(i)), s.D{1}(i), s.D{2}(i)));
	if i==length(s.k),
		xlabel('t');
	else
		set(gca,'xticklabel',[]);
	end
end

figure(s.fig1)
axes('Position',[0,(s.nsub-s.nk)/s.nsub,1,1-(s.nsub-s.nk)/s.nsub],'Visible','off')
text(0.5,0.99, sprintf('sde6.m: %d paths @ h=%g T=%g', s.ns,s.h,s.h*s.maxiter), ...
		'HorizontalAlignment','center', ...
		'VerticalAlignment','top','Units','normalized');
text(0.05,0.5,'<n(t)v^{*}(t)+n^{*}(t)v(t)>', ...
		'rotation',90,'HorizontalAlignment','center', ...
		'VerticalAlignment','top','Units','normalized');

figure(s.fig2)
axes('Position',[0,(s.nsub-s.nk)/s.nsub,1,1-(s.nsub-s.nk)/s.nsub],'Visible','off')
text(0.5,0.99, sprintf('sde6.m: %d paths @ h=%g T=%g', s.ns,s.h, s.h*s.maxiter), ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top','Units','normalized');
text(0.05,0.5,'<|v(t)|^2>', ...
    'rotation',90,'HorizontalAlignment','center', ...
    'VerticalAlignment','top','Units','normalized');

figure(s.fig3)
axes('Position',[0,(s.nsub-s.nk)/s.nsub,1,1-(s.nsub-s.nk)/s.nsub],'Visible','off')
text(0.5,0.99, sprintf('sde6.m: %d paths @ h=%g T=%g', s.ns,s.h, s.h*s.maxiter), ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top','Units','normalized');
text(0.05,0.5,'<|n(t)|^2>', ...
    'rotation',90,'HorizontalAlignment','center', ...
    'VerticalAlignment','top','Units','normalized');


%orient tall
%print -dpsc sde6


