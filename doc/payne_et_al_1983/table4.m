function table4
% function table4

%
% $Id: table4.m,v 1.1 2001/02/07 13:05:38 patricg Exp $
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
	T = 1/(2*pi/10);
	h=[0.0 T/15000 T/10000 T/5000 T/3500 T/2500 T/2000 T/1500 T/1250];

	epsilon = zeros(size(h));
	N = zeros(size(h));
	P = zeros(size(h));
	H = zeros(size(h));

	solver=zakharov('B','N',512,'maxiter',0,'display',0,'save',0,'h',h(1));
	[h(1), epsilon(1)] = zakharov_epsilon(solver);
	[h(1), N(1), P(1), H(1)] = zakharov_NPH(solver);

	for i=2:length(h),
		solver=zakharov('B','N',512,'display',0,'save',0,'h',h(i));
		[h(i), epsilon(i)] = zakharov_epsilon(solver);
		[h(i), N(i), P(i), H(i)] = zakharov_NPH(solver);
	end

	save table4 h epsilon N P H solver
end

fprintf(1,'\n');
fprintf(1,'      h      epsilon      N           P           H\n');
for i=length(h):-1:2,
	fprintf(1,'  %.6f   %.4f   %.6f   %.4f   %.5f\n', ...
		h(i), epsilon(i), N(i), P(i), H(i));
end
fprintf(1,'  %s             %.6f   %.4f   %.5f\n','initial',N(1),P(1),H(1));

fid=fopen('table4.dat','w');
fprintf(fid,'$h$ & $\\epsilon$ & $N$ & $P$ & $H$\\\\\n');
fprintf(fid,'\\hline\n');
for i=length(h):-1:2,
  fprintf(fid,'$%.6f$&$%.4f$&$%.6f$&$%.4f$&$%.5f$\\\\\n', ...
	    h(i), epsilon(i), N(i), P(i), H(i));
end
fprintf(fid,'%s&&$%.6f$&$%.4f$&$%.5f$\\\\\n','initial',N(1),P(1),H(1));
fclose(fid);

