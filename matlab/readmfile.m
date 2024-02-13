function params = readmfile(dpath,mfile,params)
% function params = readmfile(dpath,mfile,params)
% 
% Read content of `mfile' into structure `params'.
% `mfile' is located in `dpath' and can be generated
% from a C++ program using saveMatlab() global template
% function.

%
% $Id: readmfile.m,v 1.2 2011/03/26 09:20:41 patrick Exp $
%
% Copyright (c) 2010-2011 Patrick Guio <patrick.guio@gmail.com>
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

if ~exist('params'),
  params = [];
end

addpath(dpath)
if exist(mfile)==2, % m-file
  eval(mfile);
end
rmpath(dpath)

vars = whos;

for i=1:length(vars),
  if ~strcmp(vars(i).name,'ans') & ...
     ~strcmp(vars(i).name,'dpath') & ...
     ~strcmp(vars(i).name,'mfile') & ...
		 ~strcmp(vars(i).name,'params'),
    params = setfield(params,vars(i).name,eval(vars(i).name));
	end
end


