function stint
%
% stint  Approximate stochastic integrals
%
% Ito and Stratonovich integrals of exp(at) dW

%
% $Id: stint.m,v 1.2 2011/03/26 09:20:42 patrick Exp $
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

randn('state',100)                      % set the state of randn
T = 1; 
N = 500; 
dt = T/N; 
t = [0,dt:dt:1];

a = 100;

dB = sqrt(dt)*randn(1,N);               % increments
B = [0,cumsum(dB)];                         % cumulative sum

ito = exp(-a*t(end))*sum(exp(a*t(1:end-1)).*dB);
strat = exp(-a*t(end))*sum(exp(a*0.5*(t(1:end-1)+t(2:end))).*dB);

fprintf(1,'ito          = %f\n', ito);
fprintf(1,'stratonovich = %f\n', strat);


ito = B(end)-a*exp(-a*t(end))*sum(B(1:end-1).*exp(a*t(1:end-1))*dt);
strat = B(end)-a*exp(-a*t(end))*sum(0.5*(B(1:end-1)+B(2:end)).* ...
	exp(a*0.5*(t(1:end-1)+t(2:end)))*dt);

fprintf(1,'ito          = %f\n', ito);
fprintf(1,'stratonovich = %f\n', strat);
