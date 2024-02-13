function fig9
% function fig9

%
% $Id: fig9.m,v 1.2 2001/02/08 18:52:33 patricg Exp $
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

h=subplot(1,1,1);
lim=get(h,'Position'); 
left=lim(1); bottom=lim(2); width=lim(3); height=lim(4);
clf

axes('position', [left,bottom+2*height/3,width,height/3]);
plot(solver.xj, real(solver.Ej))
set(gca,'xticklabel',[]);
set(gca,'xlim',[-80 80],'ylim',[-0.2 0.2]);
ticks=get(gca,'yticklabel'); ticks(1,:)=' '; set(gca,'yticklabel',ticks);
legend('REAL',4)

axes('position', [left,bottom,width,2*height/3]);
plot(solver.xj, imag(solver.Ej))
set(gca,'xlim',[-80 80],'ylim',[-1.0 0.2]);
xlabel('x')
legend('IMAGINARY',4)

orient landscape; set(gcf,'PaperOrientation','portrait');
if exist('exportfig'),
	exportfig(gcf,'fig9');
else
	print -deps fig9
end

