function solver=efluct_init(solver)
% function solver=efluct_init(solver)

%
% $Id: efluct_init.m,v 1.8 2011/03/26 09:20:41 patrick Exp $
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

solver.nk = zeros(size(solver.k));
solver.dnk = zeros(size(solver.k));

solver.Ek = solver.Erms*exp(i*2*pi*rand(size(solver.k)));
solver.Ek(find(solver.k==0)) = solver.E0;

% Set to zero the k-modes larger than M
solver.Ek(solver.ikeq0) = 0.0;
solver.nk(solver.ikeq0) = 0.0;

% Fourier transform 
solver.Ej = IFFT(solver.Ek, solver.N, solver.M);
solver.nj = IFFT(solver.nk, solver.N, solver.M);
solver.dnj = IFFT(solver.dnk, solver.N, solver.M);

