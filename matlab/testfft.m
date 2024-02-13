function testfft
% function testfft

%
% $Id: testfft.m,v 1.3 2011/03/26 09:20:42 patrick Exp $
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

n=64;
test(n)

n=65;
test(n)


function test(n)
index=[-fix(n/2):fix(n/2)-1]';
x=ones(1,n);
X=fftshift(fftw(x,-1));
ii=find(index==0);
if X(ii)==n
	fprintf(1,'n=%d OK\n',n);
else
	fprintf(1,'n=%d not OK\n',n);
end

