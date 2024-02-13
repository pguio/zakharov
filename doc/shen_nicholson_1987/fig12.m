function fig12
% function fig12

%
% $Id: fig12.m,v 1.6 2004/04/27 07:05:51 patricg Exp $
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

addpath ../../hdf

filename = 'shen_nicholson_1987_fig12.hdf';

if exist(filename,'file');
	filename = which(filename);
	Er=sd(filename,'re-Ej');
	Ei=sd(filename,'im-Ej');
	E = sqrt(Er.^2+Ei.^2);
	E=struct(E);
	t=E.dims{1};
	Ej=E.var;
	[d,it]=min(abs(t-1));

	Ej = Ej(it,:)';
	n=sd(filename,'re-nj');
	n=struct(n);
	nj=n.var;
	xj=n.dims{2};
	nj = nj(it,:)';

	xlim=[min(xj) max(xj)];

	subplot(211),
	plot(xj,Ej)
	ylabel('|E(x)|')
	title(sprintf('Zakharov E0=10.0 T=%d',t(it)));
	set(gca,'xlim',xlim)
	ylim = [min(Ej) max(Ej)];
	set(gca,'ylim',ylim);

	Ej2 = Ej.*Ej;
	Ej2=Ej2(:);
	xj=xj(:);

	W0=integrate(xj,Ej2)/(max(xj)-min(xj));

	subplot(212),
	plot(xj,nj,'-', xj, -Ej2+W0,'--')
	xlabel('position x');
	ylabel('|n(x)|')
	set(gca,'xlim',xlim)
	set(gca,'ylim',[-0.035 0.015]);
	ylim = [min([nj;-Ej2+W0]) max([nj;-Ej2+W0])];
	set(gca,'ylim',ylim);

	%orient landscape; set(gcf,'PaperOrientation','portrait');
	orient tall
	if exist('exportfig'),
	  exportfig(gcf,'fig12');
	else
	  print -deps fig12
	end

end
