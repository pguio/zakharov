function testall
% function testall

%
% $Id: testall.m,v 1.1 2001/02/07 13:05:38 patricg Exp $
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

for i=[1:4],
  disp(['table' int2str(i)]);
	eval(['table' int2str(i)]);
end

for i=[1:5 7:12],
  disp(['fig' int2str(i)]);
  eval(['fig' int2str(i)]);
end

close all
