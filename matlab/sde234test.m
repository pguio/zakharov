function sde234test
% function sde234test

%
% $Id: sde234test.m,v 1.3 2011/03/26 09:20:42 patrick Exp $
%
% Copyright (c) 2001-2011 Patrick Guio <patrick.guio@gmail.com>
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

create = 1;
for k=[0 0.5 1.0 2.0],
	clear s
	s.k = k;
	s=sde2(s);
	orient landscape
	if create
		print -dpsc sde234test
		create = ~create;
	else
		print -dpsc -append sde234test
	end
	s=sde3(s);
	orient landscape
	if create
		print -dpsc sde234test
		create = ~create;
	else
		print -dpsc -append sde234test
	end
	s=sde4(s);
	orient landscape
	if create
		print -dpsc sde234test
		create = ~create;
	else
		print -dpsc -append sde234test
	end
end
