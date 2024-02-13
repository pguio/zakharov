function y = dEk2dp(x,f,p,dp,func)
% function y = dEk2dp(x,f,p,dp,func)

%
% $Id: dEk2dp.m,v 1.1 2001/02/08 18:52:34 patricg Exp $
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
y=[sech(p(2)*x) p(1)*(-2.0*sech(p(2)*x).^2.*tanh(p(2)*x)*p(2))];
