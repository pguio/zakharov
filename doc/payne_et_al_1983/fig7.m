function fig7
% function fig7

%
% $Id: fig7.m,v 1.2 2001/02/08 18:52:33 patricg Exp $
%
% Copyright (c) 2000 Patrick Guio <patrick@phys.uit.no>
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

close all

if exist('fig5_10.mat')
	load fig5_10
else
	solver = zakharov('C','display',0);
	save fig5_10 solver
end

subplot(1,1,1);
lim=get(gca,'Position');
left=lim(1); bottom=lim(2); width=lim(3); height=lim(4);
clf reset

d = 1e-3;

axes('position', [left,bottom+6*height/7,width,height/7]);
[mn, imn]=min(solver.times);
plot(solver.xj,abs(solver.Ejs(:,imn)))
set(gca,'xticklabel',[]);
set(gca,'xlim',[-80 80],'ylim',[0 1.2]);
ticks=get(gca,'yticklabel'); ticks(1,:)=' '; set(gca,'yticklabel',ticks);
legend(sprintf('%.1f',solver.times(imn)))

axes('position', [left,bottom+5*height/7,width,height/7]);
[mn, imn]=min(abs(solver.times-12.7));
plot(solver.xj,abs(solver.Ejs(:,imn)))
set(gca,'xticklabel',[]);
set(gca,'xlim',[-80 80],'ylim',[0 1.2]);
ticks=get(gca,'yticklabel'); ticks(1,:)=' '; set(gca,'yticklabel',ticks);
legend(sprintf('%.1f',solver.times(imn)))

axes('position', [left,bottom+3*height/7,width,2*height/7]);
[mn, imn]=min(abs(solver.times-15.9));
plot(solver.xj,abs(solver.Ejs(:,imn)))
set(gca,'xticklabel',[]);
set(gca,'xlim',[-80 80],'ylim',[0 2.4]);
ticks=get(gca,'yticklabel'); ticks(1,:)=' '; set(gca,'yticklabel',ticks);
legend(sprintf('%.1f',solver.times(imn)))

axes('position', [left,bottom+2*height/7,width,height/7]);
[mn, imn]=min(abs(solver.times-19.1));
plot(solver.xj,abs(solver.Ejs(:,imn)))
set(gca,'xticklabel',[]);
set(gca,'xlim',[-80 80],'ylim',[0 1.2]);
legend(sprintf('%.1f',solver.times(imn)))

axes('position', [left,bottom+height/7,width,height/7]);
[mn, imn]=min(abs(solver.times-25.5));
plot(solver.xj,abs(solver.Ejs(:,imn)))
set(gca,'xticklabel',[]);
set(gca,'xlim',[-80 80],'ylim',[0 1.2]);
ticks=get(gca,'yticklabel'); ticks(1,:)=' '; set(gca,'yticklabel',ticks);
legend(sprintf('%.1f',solver.times(imn)))

axes('position', [left,bottom,width,height/7]);
[mn, imn]=min(abs(solver.times-31.8));
plot(solver.xj,abs(solver.Ejs(:,imn)))
set(gca,'xlim',[-80 80],'ylim',[0 1.2]);
xlabel('x')
ticks=get(gca,'yticklabel'); ticks(1,:)=' '; set(gca,'yticklabel',ticks);
legend(sprintf('%.1f',solver.times(imn)))

orient landscape; set(gcf,'PaperOrientation','portrait');
if exist('exportfig'),
	exportfig(gcf,'fig7');
else
	print -deps fig7
end

