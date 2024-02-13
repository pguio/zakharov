function [t,N,P,H]=zakharov_NPHs(solver)
% function [t,N,P,H]=zakharov_NPHs(solver)

%
% $Id: zakharov_NPHs.m,v 1.4 2011/03/26 09:20:42 patrick Exp $
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

s = solver;

t = [solver.times solver.ctime];
N = zeros(size(t));
P = zeros(size(t));
H = zeros(size(t));

for i=1:length(t)-1,
	s.Ej = solver.Ejs(:,i);
	s.nj = solver.njs(:,i);
	s.dnj = solver.dnjs(:,i);
	s.dnk = FFT(s.dnj, s.N, s.M);
	s.Ek = FFT(s.Ej, s.N, s.M);
  	[h, N(i), P(i), H(i)] = zakharov_NPH(s);
end

[h, N(i+1), P(i+1), H(i+1)] = zakharov_NPH(solver);

