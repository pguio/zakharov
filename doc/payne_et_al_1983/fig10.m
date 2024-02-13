function fig10
% function fig10

%
% $Id: fig10.m,v 1.3 2001/09/27 07:48:39 patricg Exp $
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

subplot(111)
plot(solver.xj, abs(solver.Ej),'-',solver.xj, real(solver.nj),'--')
set(gca,'xlim',[-30 30],'ylim',[-1.5 1.0]);
xlabel('x')
legend('|E|','n',4)


orient landscape; set(gcf,'PaperOrientation','portrait');
if exist('exportfig'),
	exportfig(gcf,'fig10');
else
	print -deps fig10
end

