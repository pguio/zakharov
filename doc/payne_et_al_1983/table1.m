function table1
% function table1

%
% $Id: table1.m,v 1.1 2001/02/07 13:05:37 patricg Exp $
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

if exist('table1.mat')
  load table1
else
	sa = zakharov('A','maxiter',0,'display',0,'save',0);
	sb = zakharov('B','maxiter',0,'display',0,'save',0);
	sc = zakharov('C','maxiter',0,'display',0,'save',0);

	save table1 sa sb sc
end

fprintf(1,'\n');
fprintf(1,'      L    Emax     Emin        v           u          n0\n');
fprintf(1,'A   %5.1f  %4.1f  %.4e %10.6f %10.5f %10.6f\n', ...
	sa.L,sa.Emax, sa.Emin, sa.v, sa.u, sa.n0);
fprintf(1,'B   %5.1f  %4.1f  %.4e %10.6f %10.3f %10.5f\n', ...
	sb.L,sb.Emax, sb.Emin, sb.v, sb.u, sb.n0);
fprintf(1,'C   %5.1f  %4.1f  %.4e %10.6f %10.5f %10.7f\n', ...
	sc.L,sc.Emax(1), sc.Emin(1), sc.v(1), sc.u(1), sc.n0(1));
fprintf(1,'C   %5.1f  %4.1f  %.4e %10.6f %10.5f %10.7f\n', ...
	sc.L,sc.Emax(2), sc.Emin(2), sc.v(2), sc.u(2), sc.n0(2));

fid=fopen('table1.dat','w');
fprintf(fid,'& $L$ & $E_{max}$ & $E_{min}$ & $v$ & $u$ & $n_0$\\\\\n');
fprintf(fid,'\\hline\n');
fprintf(fid,'$A$&$%5.1f$&$%4.1f$&$%.4e$&$%10.6f$&$%10.5f$&$%10.6f$\\\\\n', ...
  sa.L,sa.Emax, sa.Emin, sa.v, sa.u, sa.n0);
fprintf(fid,'$B$&$%5.1f$&$%4.1f$&$%.4e$&$%10.6f$&$%10.3f$&$%10.5f$\\\\\n', ...
	sb.L,sb.Emax, sb.Emin, sb.v, sb.u, sb.n0);
fprintf(fid,'$C$&$%5.1f$&$%4.1f$&$%.4e$&$\\pm%10.6f$&$\\mp%10.5f$&$%10.7f$\\\\\n', ...
	sc.L,sc.Emax(1), sc.Emin(1), sc.v(1), -sc.u(1), sc.n0(1));

fclose(fid);
