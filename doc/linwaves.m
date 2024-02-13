function s=linwaves
% function s=linwaves

%
% $Id: linwaves.m,v 1.20 2023/02/23 08:54:48 patrick Exp $
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

dofig = [0 0 0 1 0 0 0 0];
dofig = [1 0 0 0 0 0 0 0];

%orient landscape; set(gcf,'PaperOrientation','portrait');
% Bounds [tight]/loose
% BoundsCode [internal]/mcode
% LockAxes [1]/0
% FontMode [scaled]/fixed
% LineMode []/scaled/fixed
opts = struct('Color','cmyk',...
              'Bounds','tight',...
              'BoundsCode','mcode',...
              'FontMode','scaled','FontSize',1.5,...
              'LineMode','scaled','LineWidth',1.5,...
              'LockAxes',1);


if dofig(1) | ~exist('linwaves1.eps','file'),
	s.seeds = 1:500;
	s.h = 1e-3;
	s.maxiter = 10000;
	s.nk = 2;
	s.k = [25 80];
	s.k = [80 80];
	s.nu = [.01 400];
	%s.nu = [1 1];
	s.phi0 = [0.002 0.002];
	s.phi0 = [1 1];
	s.D = sqrt(s.nu).*abs(s.phi0);
	%s.D = 1/sqrt(2)*sqrt(s.nu).*abs(s.phi0);
	s.nsub = 3;
	s = sde1(s);
	%subplot(3,1,1), set(gca,'ylim',[3,5])
	%subplot(3,1,2), set(gca,'ylim',[3,5])
  set(gcf,'Position',[1400 400 1200 900]);
	orient landscape, set(gcf,'PaperOrientation','portrait');
	%orient tall
	exportfig(gcf,'linwaves1',opts);
end

if dofig(2) | ~exist('linwaves2.eps','file'),
	clear s
	s.seeds = 1:300;
	s.h = 1e-3;
	s.maxiter = 1000;
	s.nk = 2;
	s.k = [25 80];
	s.nu = [4 15];
	theta = acos(-sqrt(4*s.nu.^2.*s.k.^2./(s.k.^2.*(s.k.^2+4*s.nu.^2))));
	s.phi0 = 1.5e6*ones(1,s.nk).*exp(i*theta/2);
	s.D = { sqrt(2*s.nu.*s.k.^2./(4*s.nu.^2+s.k.^2)).*abs(s.phi0), zeros(1,s.nk) };
	s.dphi0 = s.D{1}.*abs(s.k)./sqrt(2*s.nu).*exp(-i*theta/2); 
	s.nsub = 3;
	s = sde5(s);
	figure(s.fig3)
  set(gcf,'Position',[1400 400 1200 900]);
	orient landscape, set(gcf,'PaperOrientation','portrait');
	%orient tall
	exportfig(gcf,'linwaves2',opts);
end

if dofig(3) | ~exist('linwaves3.eps','file'),
	clear s
	s.seeds = 1:300;
	s.h = 1e-3;
	s.maxiter = 1000;
	s.nk = 2;
	s.k = [25 80];
	s.nu = [4 15];
	theta = acos(0)*ones(1,s.nk);
	s.phi0 = 1.5e6*ones(1,s.nk).*exp(i*theta/2);
	s.D = { zeros(1,s.nk), sqrt(2*s.nu.*s.k.^2).*abs(s.phi0) };
	s.dphi0 = s.D{2}./sqrt(2*s.nu).*exp(-i*theta/2);
	s.nsub = 3;
	s = sde6(s);
	figure(s.fig3)
  set(gcf,'Position',[1400 400 1200 900]);
	orient landscape, set(gcf,'PaperOrientation','portrait');
	%orient tall
	exportfig(gcf,'linwaves3',opts);
end

maxiter = 1000;
seeds = 1;

if dofig(4) | ~exist('linwaves4.eps','file'),
	clear s
	s.seeds=seeds;
	s.maxiter = maxiter;
	s.L_r = 25;
	s.N = 1024;
	s=linearwaves(s);
	figure(1)
	%subplot(211), set(gca,'ylim',[0 2e7]);
	%subplot(212), set(gca,'ylim',[0 0.6e-12]);
  set(gcf,'Position',[1400 400 1200 900]);
	orient landscape, set(gcf,'PaperOrientation','portrait');
	exportfig(gcf,'linwaves4',opts);
end

if dofig(5) | ~exist('linwaves5.eps','file'),
	clear s
	s.seeds=seeds;
	s.L_r = 25;
	s.N = 1024;
	s.c=2;
	s.maxiter = maxiter*s.c;
	s=linearwaves(s);
	figure(1)
%	subplot(211), set(gca,'ylim',[0 2e7]);
%	subplot(212), set(gca,'ylim',[0 0.6e-12]);
  set(gcf,'Position',[1400 400 1200 900]);
	orient landscape, set(gcf,'PaperOrientation','portrait');
	exportfig(gcf,'linwaves5',opts);
end

if dofig(6) | ~exist('linwaves6.eps','file'),
	clear s
	s.seeds=1;
	s.L_r = 25;
	s.N = 1024;
	s.c=5;
	s.maxiter = maxiter*s.c;
	s=linearwaves(s);
	figure(1)
  set(gcf,'Position',[1400 400 1200 900]);
	orient landscape, set(gcf,'PaperOrientation','portrait');
	exportfig(gcf,'linwaves6','color','cmyk');
end

if dofig(7) | ~exist('linwaves4alt.eps','file'),
	clear s
	s.seeds=1;
	s.maxiter = maxiter;
	s=linearwaves1(s);
  set(gcf,'Position',[1400 400 1200 900]);
	orient landscape, set(gcf,'PaperOrientation','portrait');
	exportfig(gcf,'linwaves4alt',opts);
end

if dofig(8) | ~exist('linwaves7.eps','file'),
	clear s
	s.seeds=1;
	s.maxiter=0;
	s.c=1;
	s.maxiter = maxiter*s.c;
	s=linearwaves(s);
	figure(1);
	subplot(211),
	plot(s.k(s.ikc),s.nui(s.ikc));
	set(gca,'xlim',[min(s.k(s.ikc)) max(s.k(s.ikc))]);
	set(gca,'xticklabel',[]);
	title('\nu_i(k)')
	subplot(212),
	plot(s.k(s.ikc),s.nue(s.ikc), s.k(s.ikc),s.k(s.ikc).^2);
	set(gca,'xlim',[min(s.k(s.ikc)) max(s.k(s.ikc))]);
	title('\nu_e(k) and k^2')
	xlabel('k')

  set(gcf,'Position',[1400 400 1200 900]);
	orient landscape, set(gcf,'PaperOrientation','portrait');
	exportfig(gcf,'linwaves7',opts);

	save linwaves7 s;
end

