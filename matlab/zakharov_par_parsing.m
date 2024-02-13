function solver=zakharov_par_parsing(solver,varargin)
% function solver=zakharov_par_parsing(solver,varargin)

%
% $Id: zakharov_par_parsing.m,v 1.6 2011/03/26 09:20:42 patrick Exp $
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


if isempty(varargin),
  return
end

argc=length(varargin);

if rem(argc,2)~=0,
	error('Extra arguments should be in pair parameter/value')
end

for i=1:2:argc,
	
	param=varargin{i};
	val=varargin{i+1};

	switch param

		case 'nui', solver.nui=val;
		case 'nue', solver.nue=val;

		case 'Emax', solver.Emax=val;
		case 'm', solver.m=val;

		case 'L', solver.L=val;
		case 'N', solver.N=val;
		case 'M', solver.M=floor(val);

		case 'D', solver.D=val;

		case 'h', solver.h=val;

		case 'maxiter', solver.maxiter=val;

		case 'display', solver.display=val;
		case 'mdisplay', solver.mdisplay=val;

		case 'save', solver.save=val;
		case 'msave', solver.msave=val;
	end

end
