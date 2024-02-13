function fig3
% function fig3


%
% $Id: fig3.m,v 1.2 2001/02/08 18:52:33 patricg Exp $
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

if exist('table2.mat')
  load table2
else
	table2
	load table2
end

if exist('table3.mat')
  load table3
else
	table3
	load table3
end

h = [h64(2:end-1) h128(2:end-1)];
epsilon = [epsilon64(2:end-1) epsilon128(2:end-1)];

hi = logspace(log10(1.0e-2),log10(3.0e-1),10);
[p,s]=polyfit(log10(h),log10(epsilon),1);
epsiloni=10.^(polyval(p,log10(hi)));

h=loglog(h64, epsilon64,'o', h128, epsilon128,'p',hi,epsiloni,'-');
set(h(:),'MarkerSize',5)
xlabel('TIME STEP h');
ylabel('ERROR \epsilon');
set(gca,'xlim',[1e-3 1e0],'ylim',[1e-3 1e0]);
legend(sprintf('M=%d (N=%d)',solver64.M, solver64.N), ...
	sprintf('M=%d (N=%d)',solver128.M, solver128.N),2);

orient landscape; set(gcf,'PaperOrientation','portrait');
if exist('exportfig'),
	exportfig(gcf,'fig3');
else
	print -deps fig3
end

