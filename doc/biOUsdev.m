function sdev = biOUsdev(t,X0,nu,k,s)
% function sdev = biOUsdev(t,X0,nu,k,s)
%
% Compute the standard deviation of the real and imaginary parts of the
% complex process \sqrt{<(X_r(t)-<X_r(t)>)^2>} and 
% \sqrt{<(X_i(t)-<X_i(t)>)^2>} of the complex Ornstein-Uhlenbeck 
% process defined by the complex Langevin equation
% X(t+dt) = X(t) -(nu + i *k^2) X(t) dt + dW1 + i*dW2
%

%
% $Id: biOUsdev.m,v 1.1 2010/06/16 14:29:54 patrick Exp $
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


var = biOUvar(t,X0,nu,k,s);
sdev = sqrt(real(var))+i*sqrt(imag(var));

