function fig9
%function fig9

%
% $Id: fig9.m,v 1.5 2004/04/27 07:05:51 patricg Exp $
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

clf 

addpath ../../hdf

filename = 'shen_nicholson_1987_fig8.hdf';

if exist(filename,'file');
  filename = which(filename);
	Ek2=sd(filename,'Ek2');
	Ek2=struct(Ek2);
	k=Ek2.dims{2};
	Ek2=Ek2.var;

	ii = find(k>0);
	k = k(ii); 
	k = k(:);
	Ek2 = Ek2(ii);
	Ek2 = Ek2(:);

	k_dk=k/(k(2)-k(1));


	p=[3.11e-2; 0.933];
	options=[[0.0; 0.0] [inf; inf]];
	[Ek2f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2]=...
		leasqr(k,Ek2,p,'Ek2',.0001,50,1,[1;1],'dEk2dp',options);

	fprintf(1,'<|Ek|^2>= %e sech^2(%f k)\n', p);
	fid=fopen('fig9.dat','w');
	fprintf(fid,'$|E_k|^2= %f \\sech^2(%f k)$', p);
	fclose(fid);

	subplot(111)
	loglog(k_dk, Ek2f, k_dk, Ek2, 'o')
	set(gca,'xlim',[1 100]);
	set(gca,'ylim',[1e-4 1e-1]);

	xlabel('wavenumber k/\Delta k');
	ylabel('|E(k)|^2');
	title(sprintf('Zakharov E0=1.0'));

	%orient landscape; set(gcf,'PaperOrientation','portrait');
	orient portrait;
	if exist('exportfig'),
	  exportfig(gcf,'fig9');
	else
	  print -deps fig9
	end

end

