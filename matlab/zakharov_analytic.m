function solver=zakharov_analytic(f,varargin)
% function solver=zakharov_analytic(f,varargin)

%
% $Id: zakharov_analytic.m,v 1.9 2011/03/26 09:20:42 patrick Exp $
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

close all

solver = feval(f);
solver = zakharov_par_parsing(solver,varargin{:});
solver = zakharov_init(solver);

if solver.nb~=1
	error('zakharov_analytic can only handle one soliton');
end
if solver.display
	plot_zakharov(solver);
end
view_zakharov(solver);

if solver.save
	solver.Ejs=solver.Ej;
	solver.njs=solver.nj;
	solver.dnjs=solver.dnj;
	solver.times=0;
end

N = solver.N;
M = solver.M;

for i=1:solver.maxiter,

	solver.citer = i;
	solver.ctime = solver.ctime + solver.h;
	[E,n,dn]=analytic_one_soliton(solver,1);
	solver.Ej=E;
	solver.nj=n;
	solver.dnj=dn;

	% Fourier transform 
	solver.nk = FFT(solver.nj, N, M);
	%solver.nk(find(solver.k==0))=0;
	solver.dnk = FFT(solver.dnj, N, M);
	%solver.dnk(find(solver.k==0))=0;
	solver.Ek = FFT(solver.Ej, N, M);
	%solver.Ek = imag(FFT(solver.Ej, N, M));

	if mod(i,100)==0,
		s=sprintf(' Iter %6.d/%6.d: integrated to %.2f/%.2f (%.2f%%)', ...
			i,solver.maxiter,solver.ctime,solver.maxtime, ...
			100*solver.ctime/solver.maxtime);
		fprintf(1,'%s%s', s, repmat(8,1,length(s)));
	end

	if solver.save & mod(i,solver.msave)==0,
		solver.Ejs=[solver.Ejs solver.Ej];
		solver.njs=[solver.njs solver.nj];
		solver.dnjs=[solver.dnjs solver.dnj];
		solver.times=[solver.times solver.ctime];
	end

	if solver.display & mod(i,solver.mdisplay)==0,
		plot_zakharov(solver);
	end

end

if solver.maxiter,
	fprintf(1,' Iter %6.d/%6.d: integrated to %.2f/%.2f (%.2f%%)\n', ...
  	i,solver.maxiter,solver.ctime,solver.maxtime, ...
		100*solver.ctime/solver.maxtime);
end

if solver.display
  plot_all_zakharov(solver);
end
