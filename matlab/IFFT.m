function f=IFFT(F,N,M)
% function f=IFFT(F,N,M)

%
% $Id: IFFT.m,v 1.7 2011/03/26 09:20:41 patrick Exp $
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

% set to zero the k-modes with index larger than M
signed_index=[-N/2:N/2-1]';
F(find(signed_index>M | signed_index<-M)) = 0.0;

F=ifftshift(F);

if exist('fftw')==3,
	f = fftw(F,1)/prod(size(F));
else
	f = ifftn(F);
end

