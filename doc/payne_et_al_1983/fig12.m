function fig11
% function fig11

%
% $Id: fig12.m,v 1.2 2001/02/08 18:52:33 patricg Exp $
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

if exist('fig11_12.mat')
	load fig11_12
else
	T = 20/(2*pi/10);

	h1 = T/1000;
	solver1=zakharov('C','h',h1,'display',0,'msave',5);
	[t1, N1, P1, H1] = zakharov_NPHs(solver1);

	h2 = T/2000;
	solver2=zakharov('C','h',h2,'display',0,'msave',10);
	[t2, N2, P2, H2] = zakharov_NPHs(solver2);

	save fig11_12 h1 t1 N1 H1 h2 t2 N2 H2
end

subplot(111)
plot(t1,(H1(1)-H1)/H1(1),'-', t2, (H2(1)-H2)/H2(1), '--')
set(gca,'xlim',[0 max([t1(:);t2(:)])])
set(gca,'ylim',[-0.012 0.002])
xlabel('TIME t')
ylabel('ERROR \epsilon_H');
legend(sprintf('%.4f',h1), sprintf('%.4f',h2),4)

orient landscape; set(gcf,'PaperOrientation','portrait');
if exist('exportfig'),
	exportfig(gcf,'fig12');
else
	print -deps fig12
end

