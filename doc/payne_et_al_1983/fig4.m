function fig4
% function fig4


%
% $Id: fig4.m,v 1.2 2001/02/08 18:52:33 patricg Exp $
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

if exist('table4.mat')
  load table4
else
	table4
	load table4
end

h = h(2:end);
epsilon = epsilon(2:end);

hi = logspace(log10(3.0e-5),log10(1.5e-3),5);
[p,s]=polyfit(log10(h),log10(epsilon),1);
epsiloni=10.^(polyval(p,log10(hi)));

h=loglog(h, epsilon,'o',hi, epsiloni,'-');
set(h(:),'MarkerSize',5)
xlabel('TIME STEP h');
ylabel('ERROR \epsilon');
set(gca,'xlim',[1e-5 1.5e-3],'ylim',[1e-3 1.0e0]);
legend(sprintf('M=%d (N=%d)',solver.M, solver.N),2);

orient landscape; set(gcf,'PaperOrientation','portrait');
if exist('exportfig'),
	exportfig(gcf,'fig4');
else
	print -deps fig4
end
