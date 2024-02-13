function [t,epsilon]=zakharov_epsilons(solver)
% function [t,epsilon]=zakharov_epsilons(solver)
%
% Calculate epsilon given by Eq.33 of payne (1983)

%
% $Id: zakharov_epsilons.m,v 1.2 2011/03/26 09:20:42 patrick Exp $
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
epsilon = zeros(size(t));

for i=1:length(t)-1,
	s.Ej = solver.Ejs(:,i);
	s.ctime = t(i);
	[h,epsilon(i)]=zakharov_epsilon(s);
end

[h, epsilon(i+1)] = zakharov_epsilon(solver);

