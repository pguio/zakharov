function [var,varargout] = biOUvar(t,X0,nu,k,s)
% function [var[,covar]] = biOUvar(t,X0,nu,k,s)
%
% Compute the variances (and optionally covariance) of the real and
% imaginary parts of the complex Ornstein-Uhlenbeck process
% <(X_r(t)-<X_r(t)>)^2> and <(X_i(t)-<X_i(t)>)^2> defined by the complex
% Langevin equation
% X(t+dt) = X(t) -(nu + i *k^2) X(t) dt + dW1 + i*dW2
%

%
% $Id: biOUvar.m,v 1.6 2010/06/17 17:26:22 patrick Exp $
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


x0 = real(X0);
y0 = imag(X0);

m = biOUmean(t,X0,nu,k);

% variance of real part
vr = 1/2*exp(-2*nu*t).*(x0^2+y0^2+...
                        (x0^2-y0^2)*cos(sqrt(2)*k^2*t) + ...
                        sqrt(2)*x0*y0*sin(sqrt(2)*k^2*t)) + ...
    s^2/(2*nu)*(1-exp(-2*nu*t)) - real(m).^2;

% variance of imaginary part
vi = 1/2*exp(-2*nu*t).*(x0^2+y0^2-...
                        (x0^2-y0^2)*cos(sqrt(2)*k^2*t) + ...
                        sqrt(2)*x0*y0*sin(sqrt(2)*k^2*t)) + ...
    s^2/(2*nu)*(1-exp(-2*nu*t)) - imag(m).^2;

var = vr + i*vi;

% covariance between real and imaginary parts
if nargout>1,
vri = 1/2*exp(-2*nu*t).*(...%sqrt(2)*(x0^2-y0^2)*sin(sqrt(2)*k^2*t) + ...
                         2*x0*y0*cos(sqrt(2)*k^2*t)) - real(m).*imag(m);
varargout(1) = {vri};
end

