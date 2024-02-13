function linwaves4096(filename,seeds)
% function linwaves4096(filename,seeds)

%
% $Id: linwaves4096.m,v 1.7 2011/03/26 09:20:41 patrick Exp $
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


if nargin==0,
	d=sd('hdf/linwaves4096/linwaves4096','dkf2',1);
	i=sd('hdf/linwaves4096/linwaves4096','nkf2',1);
	u=sd('hdf/linwaves4096/linwaves4096','ukf2',1);
else
	d=sd(filename,'dkf2',[],seeds);
	i=sd(filename,'nkf2',[],seeds);
	u=sd(filename,'ukf2',[],seeds);
end



imagesc(d', [2 3 1])  
axis ij
imagesc(i', [2 3 2])  
imagesc(u', [2 3 3])  

[KD,D] = specint(d);
[KI,I] = specint(i);
[KU,U] = specint(u);

subplot(223), plot(KI, I)
subplot(224), plot(abs(KD), D, KU, U);


function [k,i]=specint(s)
s = struct(s);
S = s.var;
k = s.dims{1};
[F,K] = meshgrid(s.dims{2}, s.dims{1});
i = integrate(F', S');

