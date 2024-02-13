function view_zakharov(solver)
% function view_zakharov(solver)

%
% $Id: view_zakharov.m,v 1.9 2011/03/26 09:20:42 patrick Exp $
%
% Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
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
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

fprintf(1,'*************************************\n');
fprintf(1,'1-D PERIODIC ZAKHAROV EQUATION SOLVER\n');
fprintf(1,'*************************************\n');

fprintf(1,'SPATIAL PERIOD L=%.2f (%.2f m)\n', solver.L, solver.L*solver.chi);
fprintf(1,'DISTANCE D=%.2f (%.2f m)\n', solver.D, solver.D*solver.chi);
fprintf(1,'GRID SIZE CONFIGURATION SPACE N=%d\n', solver.N);
fprintf(1,'TRUNCATION AT MODE M=%d\n', solver.M);

fprintf(1,'TIME INCREMENT h=%.2e (%.2e s)\n', solver.h, solver.h*solver.tau);

fprintf(1,'MAXITER=%d\n', solver.maxiter);

if isfield(solver,'solitons')
	for k=1:solver.nb,
		fprintf(1,'SOLITON #%d\n', k);
		fprintf(1,'\t=> EMAX=%f EMIN=%e\n', solver.Emax(k), solver.Emin(k));
		fprintf(1,'\t=> M=%d V=%f U=%f N0=%f\n', ...
			solver.m(k), solver.v(k), solver.u(k), solver.n0(k));
		fprintf(1,'\t=> Q=%e\n', solver.q(k));
	end
end

fprintf(1,'*************************************\n');
