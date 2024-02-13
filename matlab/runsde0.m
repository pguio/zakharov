function runsde0
% function runsde0

%
% $Id: runsde0.m,v 1.6 2011/03/26 09:20:41 patrick Exp $
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

s.c=1/50; 
s.seeds=[1:1000];
s=sde0(s); 
orient landscape
print -dpsc sde0
for c=[1/5 1/2.5 1 2 4],
	clear s
	s.seeds=[1:1000];
	s.c=c;
	s=sde0(s);
	orient landscape
	print -dpsc -append sde0
end

