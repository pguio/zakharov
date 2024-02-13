function s = nyesub(s)
% function s = nyesub(s)
%
% Langmuir wavenumber damping operator from linear theory 

%
% $Id: nyesub.m,v 1.3 2011/03/26 09:20:41 patrick Exp $
%
% Copyright (c) 2001-2011 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

eta = 1+3/s.t_ratio;

s1  = sqrt(pi/8)*1.5^4*exp(-1.5)*(s.mass_ratio/eta)^2.5;
s2  = 9*s.mass_ratio/(8*eta);

if 0
  kbreak = s.nbreak*s.pi2ndx;
	k2break = kbreak^2;
	gammabreak = (s1/(k2break*kbreak))*exp(-s2/kbreak);

	beta = (2*s2/k2break-3)*gammabreak/(2*k2break);
	alpha = gammabreak-beta*k2break;

	if s.nbreak>s.n2-1 | beta<0,
	  error('Cannot assign this Landau damping extension');
	end
end

nye = zeros(size(1:s.n));
m=1:s.n2-1;
k2 = s.pi2ndx2*m.^2;
nye(m) = (s1./(k2.*s.pi2ndx.*m)).*exp(-s2./k2) + s.ecoll/2;
nye(s.n+1-m) = nye(m);

if 0
  m=s.nbreak:s.n2-1;
  k2 = s.pi2ndx2*m.^2;
  nye(m) = alpha + beta*k2 + ecoll/2;
  nye(s.n+1-m) = nye(m);
end

% k=0 mode
nye(1) = s.ecoll/2;

% k=-(n/2)dk mode
m=s.n2+1;
k2 = s.pi2ndx2*m.^2;
if 0
  nye(m) = alpha + beta*k2 + ecoll/2;
else
  nye(m) = (s1/(k2*s.pi2ndx*s.n2))*exp(-s2/k2) + s.ecoll/2;
end


s.nye = nye(:);
