function dy=derivate(x,y)
% function dy=derivate(x,y)

%
% $Id: derivate.m,v 1.5 2011/03/26 09:20:41 patrick Exp $
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

if 0 & exist('derivspl')==3,

	dy = derivspl(x(:), y(:), x(:));

%	dy1 = zeros(size(y(:)));
%	dy1(1) = (y(end)-y(2)) / 2.0 / (x(2)-x(1));
%	dy1(2:end-1) = (y(3:end)-y(1:end-2))./ (x(3:end)-x(1:end-2));
%	dy1(end) = (y(end-1)-y(1)) / 2.0 / (x(end)-x(end-1));
%	plot(x,dy,x,dy1)
%	pause

else

	dy = zeros(size(y(:)));
	dy(1) = (y(end)-y(2)) / 2.0 / (x(2)-x(1));
	dy(2:end-1) = (y(3:end)-y(1:end-2))./ (x(3:end)-x(1:end-2));
	dy(end) = (y(end-1)-y(1)) / 2.0 / (x(end)-x(end-1));

end
