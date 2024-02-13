function fig1
% function fig1

%
% $Id: fig1.m,v 1.2 2001/02/08 18:52:33 patricg Exp $
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

if exist('fig1.mat')
 	load fig1
else
	solver=zakharov('A','N',256,'maxiter',0,'display',0,'save',0);
	x=solver.xj;
	E1=solver.Ej;
	n=solver.nj;

	solver=zakharov('B','N',256,'maxiter',0,'display',0,'save',0);
	E10=solver.Ej;

	save fig1 x E1 n E10
end

subplot(111)
plot(x,abs(E1),'-',x,abs(E10),'-',x,n,'--');
set(gca,'xlim',[-10 10],'ylim',[-2 12]);
xlabel('x');
legend('|E|','|E|','n')

orient landscape; set(gcf,'PaperOrientation','portrait');
if exist('exportfig'),
	exportfig(gcf,'fig1');
else
	print -deps fig1
end

