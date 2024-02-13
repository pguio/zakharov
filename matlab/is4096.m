function is4096
% function is4096

%
% $Id: is4096.m,v 1.4 2011/03/26 09:20:41 patrick Exp $
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

d=sd('hdf/is4096i','dkf2');
i=sd('hdf/is4096i','nkf2');
u=sd('hdf/is4096i','ukf2');

imagesc(d', [2 3 1])  
axis ij
imagesc(i', [2 3 2])  
imagesc(u', [2 3 3])  

[KD,D] = integrate(d);
[KI,I] = integrate(i);
[KU,U] = integrate(u);

subplot(223), plot(KI, I)
subplot(224), plot(abs(KD), D, KU, U);

function [k,i]=integrate(s)
s = struct(s);
S = s.var;
k = s.dims{1};
[F,K] = meshgrid(s.dims{2}, s.dims{1});
i = integrate(F', S');

