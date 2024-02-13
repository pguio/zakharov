function checkw
% function checkw
%
% Check mean and variance of complex W

%
% $Id: checkw.m,v 1.2 2011/03/26 09:20:41 patrick Exp $
%
% Copyright (c) 2007 Patrick Guio <patrick.guio@gmail.com>
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

rand('twister',1);

h=1;

for j=1:30,
	u1 = rand(1,50000);
	u2 = rand(1,50000);
	w = sqrt(h)*sqrt(-2*log(u1)).*exp(i*2*pi*u2);
	fprintf(1,'mean=%+5.3f,%+5.3f var=%5.3f, %5.3f\n', ...
		mean(real(w)), mean(imag(w)), var(real(w)), var(imag(w)));
end
