function table2
% function table2

%
% $Id: table2.m,v 1.1 2001/02/07 13:05:38 patricg Exp $
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
	T = 20/(2*pi/10);
	h64=[0.0 T/3200 T/2800 T/2400 T/2000 T/1600 T/1400 T/1200 T/1000 T/800 ...
		T/700 T/600 T/500 T/400 T/380 T/320 T/280 T/240 T/220 T/200];

	epsilon64 = zeros(size(h64));
	N64 = zeros(size(h64));
	P64 = zeros(size(h64));
	H64 = zeros(size(h64));

	solver64=zakharov('A','N',64,'maxiter',0,'display',0,'save',0,'h',h64(1));
	[h64(1), epsilon64(1)] = zakharov_epsilon(solver64);
	[h64(1), N64(1), P64(1), H64(1)] = zakharov_NPH(solver64);

	for i=2:length(h64),
		solver64=zakharov('A','N',64,'display',0,'save',0,'h',h64(i));
		[h64(i), epsilon64(i)] = zakharov_epsilon(solver64);
		[h64(i), N64(i), P64(i), H64(i)] = zakharov_NPH(solver64);
	end

	save table2 h64 epsilon64 N64 P64 H64 solver64

end

fprintf(1,'\n');
fprintf(1,'     h      epsilon      N           P           H\n');
for i=length(h64):-1:2,
	fprintf(1,'  %.5f   %.4f   %.7f   %.7f   %.7f\n', ...
		h64(i), epsilon64(i), N64(i), P64(i), H64(i));
end
fprintf(1,'  %s            %.7f   %.7f   %.7f\n','initial', ...
	N64(1), P64(1), H64(1));

fid=fopen('table2.dat','w');
fprintf(fid,'$h$ & $\\epsilon$ & $N$ & $P$ & $H$\\\\\n');
fprintf(fid,'\\hline\n');
for i=length(h64):-1:2,
	fprintf(fid,'$%.5f$&$%.4f$&$%.7f$&$%.7f$&$%.7f$\\\\\n', ...
		h64(i), epsilon64(i), N64(i), P64(i), H64(i));
end
fprintf(fid,'%s&&$%.7f$&$%.7f$&$%.7f$\\\\\n','initial', ...
	N64(1), P64(1), H64(1));
fclose(fid);

