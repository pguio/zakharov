function table3
% function table3

%
% $Id: table3.m,v 1.1 2001/02/07 13:05:38 patricg Exp $
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

if exist('table3.mat')
	load table3
else
	T = 20/(2*pi/10);
	h128=[0.0 T/1600 T/800 T/400 T/200];

	epsilon128 = zeros(size(h128));
	N128 = zeros(size(h128));
	P128 = zeros(size(h128));
	H128 = zeros(size(h128));

	solver128=zakharov('A','N',128,'maxiter',0,'display',0,'save',0,'h',h128(1));
	[h128(1), epsilon128(1)] = zakharov_epsilon(solver128);
	[h128(1), N128(1), P128(1), H128(1)] = zakharov_NPH(solver128);

	for i=2:length(h128),
		solver128=zakharov('A','N',128,'display',0,'save',0,'h',h128(i));
		[h128(i), epsilon128(i)] = zakharov_epsilon(solver128);
		[h128(i), N128(i), P128(i), H128(i)] = zakharov_NPH(solver128);
	end

	save table3 h128 epsilon128 N128 P128 H128 solver128
end

fprintf(1,'\n');
fprintf(1,'     h      epsilon      N           P           H\n');
for i=length(h128):-1:2,
	fprintf(1,'  %.5f   %.4f   %.7f   %.7f   %.7f\n', ...
		h128(i), epsilon128(i), N128(i), P128(i), H128(i));
end
fprintf(1,'  %s            %.7f   %.7f   %.7f\n','initial', ...
	N128(1), P128(1), H128(1));

fid=fopen('table3.dat','w');
fprintf(fid,'$h$ & $\\epsilon$ & $N$ & $P$ & $H$\\\\\n');
fprintf(fid,'\\hline\n');
for i=length(h128):-1:2,
	fprintf(fid,'$%.5f$&$%.4f$&$%.7f$&$%.7f$&$%.7f$\\\\\n', ...
		h128(i), epsilon128(i), N128(i), P128(i), H128(i));
end
fprintf(fid,'%s&&$%.7f$&$%.7f$&$%.7f$\\\\\n','initial', ...
	N128(1), P128(1), H128(1));
fclose(fid);

