function S=integrate(x,y)
% function S=integrate(x,y)

%
% $Id: integrate.m,v 1.5 2011/03/26 09:20:41 patrick Exp $
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

if exist('d01gaf')==3,

	S = d01gaf(x, y);

%	S1 = sum((x(2:end)-x(1:end-1))/2.0 .* (y(1:end-1) + y(2:end)));
%	fprintf(1,'S=%e S1=%e S-S1= %e\n', S, S1, S-S1);

else

	S = sum((x(2:end,:)-x(1:end-1,:))/2.0 .* (y(1:end-1,:) + y(2:end,:)));
end


