function s = sde1(s)
% function s = sde1(s)

%
% $Id: sde1.m,v 1.30 2011/03/26 09:20:42 patrick Exp $
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
	s.h = 1e-3;
end
if nargin == 0 | (nargin == 1 & ~isfield(s,'maxiter')),
	s.maxiter = 10000;
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
if nargin == 0 | (nargin == 1 & ~isfield(s,'D')),
	s.D = sqrt(s.nu).*abs(s.phi0);
end
if nargin == 0 | (nargin == 1 & ~isfield(s,'method')),
	s.method = 'exact';
%	s.method = 'euler-maruyama';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,' SDE for Langmuir waves\n');
fprintf(1,' dE = -(nu + i k^2) E dt + D dW\n\n');
fprintf(1,' nu     = '); fprintf(1,'%10.4g ', s.nu); fprintf(1,'\n');
fprintf(1,' k      = '); fprintf(1,'%10.4g ', s.k); fprintf(1,'\n');
fprintf(1,' D      = '); fprintf(1,'%10.4g ', s.D); fprintf(1,'\n');
fprintf(1,' Re E0  = '); fprintf(1,'%10.4g ', real(s.phi0)); fprintf(1,'\n');
fprintf(1,' Im E0  = '); fprintf(1,'%10.4g ', imag(s.phi0)); fprintf(1,'\n');
fprintf(1,' method = '); fprintf(1,'%s ', s.method); fprintf(1,'\n\n');

s.ns = length(s.seeds);
s.phi2s = zeros(s.maxiter, s.nk, s.ns);
s.t = [0:s.maxiter-1]'*s.h;
moduloDisplay = fix(s.ns/10);
fprintf(1,' seed = ');
for seed = s.seeds,
	s.seed = seed;
	s = solve_sde1(s);
	s.phi2s(:,:,seed) = s.phi2;
	if ~rem(seed,moduloDisplay),
		fprintf(1,'%4d/%4d%s', seed, s.ns, char(8*ones(1,9)));
	end
end
fprintf(1,'\n\n');

s = analytic_init(s);
s.phi2stats = zeros(s.maxiter,s.nk, 3);
if 0, % log-normal statistics
data = squeeze(s.phi2s(end,1,:)); % last time and first wavenumber
p = lognfit(data); % fit log-normal parameters \mu=p(1),\sigma=p(2)
fprintf(1,'log-normal fit params mu=%f sigma=%f\n', p);
fprintf(1,'log-normal m %f, arithm m %f\n',exp(p(1)+p(2)^2/2), mean(data));
fprintf(1,'log-normal s %f, arithm s %f\n',...
        sqrt((exp(p(2)^2)-1)*exp(2*p(1)+p(2)^2)),std(data));
hist(data), pause
end
s.phi2stats(:,:,1) = mean(s.phi2s,3);
onestd = std(s.phi2s,0,3);
s.phi2stats(:,:,2) = s.phi2stats(:,:,1) + onestd;
s.phi2stats(:,:,3) = s.phi2stats(:,:,1) - onestd;
s.phi2statslgd = {'m','m+s','m-s'};

plotsde(s)

function errorfill(x,y,dy)

hold on
%h=fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.85 .85 .85]);
%h=fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.8  .9  .8 ]);
h=fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.5 0.9 0.5]);
%set(h,'EdgeColor','none');
set(h,'FaceAlpha',0.5,'EdgeAlpha',0.5);
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = solve_sde1(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise \Delta W(t,k) 
rand('twister', s.seed);
u1 = rand(s.maxiter,s.nk);
u2 = rand(s.maxiter,s.nk);

switch lower(s.method),
case 'euler-maruyama',
  % Euler-Maruyama formulation
  s.wh = sqrt(s.h)*sqrt(-2*log(u1)).*exp(i*2*pi*u2);
  s.source = repmat(s.D,s.maxiter,1).*s.wh;
case 'exact',
  % Exact formulation (Gillespie, 1996)
  s.wh = sqrt(-2*log(u1)).*exp(i*2*pi*u2);
  s.source = repmat(sqrt(s.D.^2./(2*s.nu).*(1-exp(-2*s.nu*s.h))),s.maxiter,1).*s.wh;
otherwise
  error('Unknown method.')
end

% allocate phi(t,k)
s.phi = zeros(s.maxiter, s.nk);

if 0
	% trick to use filter to perform the recursion
  % phi(0) = phi0
	% phi(n+1) = phi(n)*a +  b(n) 
	b = [s.phi0; s.source(2:end,:)];
  a = exp(-(i*s.k.^2+s.nu)*s.h);
	for k=1:s.nk,
		s.phi(:,k) = filter([1],[1; -a(k)], b(:,k));
	end
	phi1 = s.phi;
end
if 1
	% exactly the same as above but without the trick
	s.phi(1,:) = s.phi0;
  a = exp(-(i*s.k.^2+s.nu)*s.h);
	for t=2:s.maxiter,
		s.phi(t,:) = s.phi(t-1,:).*a + s.source(t,:);
	end
	phi2 = s.phi;
%	if ~isfield(s,'phi1'), find(s.phi1~=s.phi); end
%	fprintf(1,'%g ', mean(abs(phi1).^2-abs(phi2).^2)); fprintf(1,'\n');
end
if 0
	% with the approximation of the exponential
	s.phi(1,:) = s.phi0;
	a = -(s.nu+i*s.k.^2)*s.h;
	for t=2:s.maxiter,
		s.phi(t,:) = s.phi(t-1,:)+ s.phi(t-1,:).*a + s.source(t,:);
	end
	phi3 = s.phi;
%	if ~isfield(s,'phi1'), find(s.phi3~=s.phi); end
%	fprintf(1,'%g ', mean(abs(phi1).^2-abs(phi3).^2)); fprintf(1,'\n');
end
%plot(s.t,abs(phi1).^2-abs(phi2).^2,'r',s.t, abs(phi1).^2-abs(phi3).^2,'g'); pause;

if 0
	% integrating the source noise before the deterministic integration
	s.phi(1,:) = s.phi0;
  a = exp(-(i*s.k.^2+s.nu)*s.h);
	for t=2:s.maxiter,
		s.phi(t,:) = (s.phi(t-1,:) + s.source(t,:)).*a;
	end
	phi4 = s.phi;
%	if ~isfield(s,'phi1'), find(s.phi1~=s.phi); end
%	fprintf(1,'%g ', mean(abs(phi1).^2-abs(phi4).^2)); fprintf(1,'\n');
end
%plot(s.t,abs(phi1).^2-abs(phi4).^2,'r'); pause;

s.phi2 = s.phi.*conj(s.phi);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = analytic_init(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t      = repmat(s.t, 1, s.nk);
nu     = repmat(s.nu, s.maxiter, 1);
phi02  = repmat(s.phi0.^2, s.maxiter, 1);
D2onu  = repmat(s.D.^2./s.nu, s.maxiter, 1);

s.phi2a = exp(-2*nu.*t).*phi02 + D2onu.*(1-exp(-2*nu.*t));


%%%%%%%%%%%%%%%%%%%
function plotsde(s)
%%%%%%%%%%%%%%%%%%%

if ~isfield(s,'nsub'),
	s.nsub = s.nk;
end

c=[0 0 1]; g=[0 .5 0]; r=[1 0 0];
for i=1:s.nk,
	subplot(s.nsub,1,i), 
	%istats = 1:6; 
	istats = 1:3; 
	%h=plot(s.t,squeeze(s.phi2stats(:,i,istats)),s.t,s.phi2a(:,i));
	%set(h(1:3),'color',c);
	%set(h(4:6),'color',g,'linestyle','--');
	%set(h(7),'color',r,'linewidth',.6);
	%set(h(4),'color',r,'linewidth',.6);
	h=plot(s.t,squeeze(s.phi2stats(:,i,1)),s.t,s.phi2a(:,i));
  %errorfill(s.t,squeeze(s.phi2stats(:,i,1)),squeeze(s.phi2stats(:,i,2)))
	set(gca,'xlim',[min(s.t) max(s.t)]);
	set(gca,'YAxisLocation','right')
	title(sprintf('\\nu=%.2g, k=%.2g, E_0=%.2g, D_E=%.2g, \\sigma_E=%.2g', ...
		s.nu(i), s.k(i), s.phi0(i), s.D(i), s.D(i)^2/s.nu(i)));
	if i==s.nk,
		xlabel('t');
		legend({'arith mean','analytic'},'location','best')
		%legend({'arith mean','analytic','1-\sigma'},'location','best')
	else
		set(gca,'xticklabel',[]);
	end
end
axes('Position',[0,(s.nsub-s.nk)/s.nsub,1,1-(s.nsub-s.nk)/s.nsub],'Visible','off')
text(0.5,1.0,sprintf('sde1.m: %d paths @ h=%g T=%g', s.ns, s.h, s.h*s.maxiter), ...
	  'HorizontalAlignment','center', ...
		'VerticalAlignment','top','Units','normalized'); 
text(0.05,0.5,'<|E(t)|^2>', ...
	'rotation',90,'HorizontalAlignment','center', ...
	'VerticalAlignment','top','Units','normalized');
%orient tall
%print -dpsc sde1
