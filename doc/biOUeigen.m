function biOUeigen(nu,k)
% function biOUeigen(nu,k)

%
% $Id: biOUeigen.m,v 1.1 2010/06/16 09:08:48 patrick Exp $
%
% Copyright (c) 2010
% Patrick Guio <p.guio@ucl.ac.uk>
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

A = [-2*nu, 0, k^2;...
     0, -2*nu, -k^2; ...
		 -k^2, k^2, -2*nu];

[V,D] = eig(A);
diag(D)
[-2*nu+i*sqrt(2)*k^2,-2*nu-i*sqrt(2)*k^2,-2*nu].'

V
V(:,1)./[1;-1;i*sqrt(2)]
V(:,2)./[1;-1;-i*sqrt(2)]
V(:,3)./[i;i;0]

%[[1;-1;i*sqrt(2)],[1;-1;-i*sqrt(2)],[1;1;0]]/sqrt(2)

