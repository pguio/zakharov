function checkeigen(nu,k)
% function checkeigen(nu,k)
%
% Check the eigen values and eigen vectors of the ion-acoustic waves formulation
% described in the linearwaves document
%

%
% $Id: checkeigen.m,v 1.2 2011/03/26 09:20:41 patrick Exp $
%
% Copyright (c) 2007-2011 Patrick Guio <patrick.guio@gmail.com>
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
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

A=[0, 0, 1; 0, -4*nu, -k^2; -2*k^2, 2, -2*nu];

De = [-2*nu+2*sqrt(nu^2-k^2),0,0;0,-2*nu,0;0,0,-2*nu-2*sqrt(nu^2-k^2)];

g = [[1, k^2-2*sqrt(nu^2-k^2)*(nu-sqrt(nu^2-k^2)),-2*(nu-sqrt(nu^2-k^2))]; ...
		[1, k^2, -2*nu]; ...
		[1, k^2+2*sqrt(nu^2-k^2)*(nu+sqrt(nu^2-k^2)),-2*(nu+sqrt(nu^2-k^2))]]';

Ve=[g(:,1)/norm(g(:,1)) g(:,2)/norm(g(:,2)) g(:,3)/norm(g(:,3))];

[Vc,Dc]=eig(A);


[Ve Vc]

[De Dc]
