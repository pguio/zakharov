function gamma=growth_rate_nse(k,E0,sgn)
% function gamma=growth_rate_nse(k,E0,sgn)

%
% $Id: growth_rate_nse.m,v 1.1 2001/02/07 16:08:14 patricg Exp $
%
% Copyright (c) 2001 Patrick Guio <patrick.guio@fys.uio.no>
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

% Eq. 13
gamma = k.*sqrt(2*E0^2-k.^2);
gamma = real(gamma);

if exist('sgn'),
	gamma = sgn*gamma;
end
