function solver=thermal_fluctuations_init(solver)
% function solver=thermal_fluctuations_init(solver)

%
% $Id: thermal_fluctuations_init.m,v 1.19 2011/03/26 09:20:42 patrick Exp $
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

si=SIconstants;

alpha=solver.k*solver.lambdae/solver.chi;
alpha2=alpha.^2;
alpha4=alpha2.^2;
Te_Ti=solver.Te/solver.Ti;

% from Sanbonmatsu et al., Quantitative comparison of reduced-description
% particle-in-cell and quasilinear-Zakharov models for parametrically
% excited Langmuir turbulence, Phys. Plasma, 7, 2824-2841, 2000
solver.delta_n=sqrt(solver.ne*(1+alpha2)./(1+alpha2+0.5*Te_Ti))* ...
	sqrt(solver.nu);

solver.delta_n=zeros(size(solver.delta_n));

if 0
solver.delta_E=sqrt(si.kb*solver.Te/si.eps0)* ...
	sqrt(1./(1+2.0*alpha2) + (1+alpha2).*alpha2./2./ ...
	(1+alpha2+0.5*Te_Ti)./(alpha2+0.5).^2)*sqrt(solver.epsilon);
else
solver.delta_E = repmat(solver.Erms,size(solver.k));
solver.delta_E(find(solver.k==0)) = solver.E0;
end


i0=find(solver.k==0);
solver.delta_n(i0)=0.0;
solver.delta_E(i0)=0.0;

if 0
subplot(211), plot(solver.k,solver.delta_n)
subplot(212), plot(solver.k,solver.delta_E)
keyboard
pause
end

if 1
subplot(211)
x=solver.k*solver.lambdae/solver.chi;
xlim=[min(x) max(x)];
plot(x, solver.delta_n*solver.nu/solver.ne);
set(gca,'xlim',xlim);
title('\deltan(k)/n_e')
subplot(212)
plot(x, solver.delta_E*solver.epsilon/sqrt(solver.ne*si.kb*solver.Te/si.eps0));
set(gca,'xlim',xlim);
title('\deltaE(k)/(n_eT_e/\epsilon_0)^{1/2}')
xlabel('k\lambda_e');
%pause
end
