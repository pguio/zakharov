function solver=C
% function solver=C

%
% $Id: C.m,v 1.8 2011/03/26 09:20:41 patrick Exp $
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
solver.L = 160.0;
% spatial grid size
solver.N = 512;
% time increment
solver.h=1e-3;
% distance D to move
solver.D=20;
% damping
solver.nui=0;
solver.nue=0;
% number of solitons 
solver.nb=2;

solver.type='solitons';
% soliton parameters
solver.Emax = [1.0 1.0];
solver.Emin = [1.0535e-31 1.0535e-31];
% paramter m eq. 31
solver.m=[8 -8];
% initial offset position 
solver.x0=[10 -10];

solver.display=1;
solver.mdisplay=100;
solver.displaytype='x';

solver.save=1;
solver.msave=100;

