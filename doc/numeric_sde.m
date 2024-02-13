function numeric_sde
% function numeric_sde

%
% $Id: numeric_sde.m,v 1.3 2010/06/12 22:21:52 patrick Exp $
%
% Copyright (c) 2001 Patrick Guio <patrick.guio@fys.uio.no>
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

seeds   = 1:300;
h       = 1e-2;
maxiter = 350;
k       = 1;
phi0    = 1;

nu      = 0.01;
method  = 'euler-maruyama';
method  = 'exact';
s = initsde(seeds,h,maxiter,nu,k,phi0,method);
s = sde1(s);
plotsde(s)
exportfig(gcf,'numeric1','color','cmyk');

clear s
nu      = 1;
method  = 'euler-maruyama';
s = initsde(seeds,h,maxiter,nu,k,phi0,method);
s = sde1(s);
plotsde(s)
exportfig(gcf,'numeric2','color','cmyk');

clear s
nu      = 10;
method  = 'euler-maruyama';
method  = 'exact';
s = initsde(seeds,h,maxiter,nu,k,phi0,method);
s = sde1(s);
plotsde(s)
exportfig(gcf,'numeric3','color','cmyk');

clear s
h       = 1e-5;
maxiter = 35000;
s.nu    = 10;
method  = 'euler-maruyama';
method  = 'exact';
s = initsde(seeds,h,maxiter,nu,k,phi0,method);
s = sde1(s);
plotsde(s)
exportfig(gcf,'numeric4','color','cmyk');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=initsde(seeds,h,maxiter,nu,k,phi0,method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s.seeds = seeds;
s.h = h;
s.maxiter = maxiter;
s.nu = nu;
s.k = k;
s.nk = length(k);
s.phi0 = phi0;
s.sigma = sqrt(s.nu).*s.phi0;
s.method = method;


%%%%%%%%%%%%%%%%%%%
function plotsde(s)
%%%%%%%%%%%%%%%%%%%

orient landscape, set(gcf,'PaperOrientation','portrait');
orient tall
close all

h=subplot(211);
ax = get(h,'Position');
set(h,'Position',[ax(1), ax(2)-.05, ax(3), ax(4)]);
plot(s.t,squeeze(s.phi2stats(:,1,:)),s.t,s.phi2a(:,1),'linewidth',1)
set(gca,'xlim',[min(s.t) max(s.t)]);
set(gca,'ylim',[0 2]);
title(sprintf('\\nu=%g, k=%g, E(0)=%g, D=%g, \\sigma=%g', ...
		s.nu, s.k, s.phi0, s.sigma, s.sigma^2/s.nu));
xlabel('t');
ylabel('<|E(t)|^2>')
axes('Position',[0 0 1 1],'Visible','off')
text(0.5,1,sprintf('%d paths @ h=%g', length(s.seeds), s.h), ...
	  'HorizontalAlignment','center', ...
		'VerticalAlignment','top','Units','normalized'); 

