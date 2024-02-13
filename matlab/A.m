function solver=A
% function solver=A

%
% $Id: A.m,v 1.7 2011/03/26 09:20:41 patrick Exp $
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

% spatial period
solver.L=20;
% spatial grid size
solver.N=64;
% time increment
solver.h=5e-3;
% distance D to move
solver.D=20;
% damping
solver.nui=0;
solver.nue=0;

solver.type='solitons';
% number of solitons 
solver.nb=1;

% soliton parameter
solver.Emax = 1.0;
% paramter m eq. 31
solver.m=1;
% initial offset position
solver.x0=0;

solver.display=1;
solver.mdisplay=100;
solver.displaytype='x';

solver.save=1;
solver.msave=100;
