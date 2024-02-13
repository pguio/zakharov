function solver=diffusion_init(solver)
% function solver=diffusion_init(solver)

%
% $Id: diffusion_init.m,v 1.22 2011/03/26 09:20:41 patrick Exp $
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

solver.vj=zeros(size(solver.k));
solver.Femj=zeros(size(solver.k));
solver.Febj=zeros(size(solver.k));
solver.Fej=zeros(size(solver.k));
solver.dFej=zeros(size(solver.k));
solver.gammaj=zeros(size(solver.k));
solver.nuej=zeros(size(solver.k));

ii=find(solver.k~=0.0);

k=solver.k/solver.chi;
tau=solver.tau;
solver.vj(ii)=solver.wpe./k(ii);

% Background
wpe=solver.wpe;
ne=solver.ne;
ve=solver.ve;
ve2=solver.ve2;
ke=1/solver.lambdae;

solver.Femj(ii)=ne/sqrt(2*pi)/ve*exp(-(ke./k(ii)).^2/2);

% Beam
if isfield(solver,'nb_n0'),
	nb=solver.ne*solver.nb_n0;
	vb=solver.vb;
	dvb=abs(solver.vb)*solver.dvb_vb;

	solver.Febj(ii)=nb/sqrt(pi)/dvb*exp(-((wpe-k(ii)*vb)./(k(ii)*dvb)).^2);
end
%
solver.Fej=solver.Femj+solver.Febj;
solver.Fej0=solver.Femj+solver.Febj;
% dFe/dv
% analytic
%dFemj(ii)=solver.Femj(ii).*(-wpe./k(ii)/ve^2);
%dFebj(ii)=solver.Febj(ii).*(-2.*(wpe-k(ii)*vb)./k(ii)/dvb^2);
% numeric
%solver.dFej=derivate(solver.vj,solver.Fej);
solver.dk=k(2)-k(1);
solver.dk_dv=-k.^2/wpe;

solver.dFej=derivate(solver.dk,solver.Fej).*solver.dk_dv;

%subplot(111), plot(k,solver.dFej,k,solver.dFej1),pause

solver.gammaj(ii)=pi/2/ne*wpe^3./k(ii).^2.*sign(k(ii)).*solver.dFej(ii)*tau;

solver.nuej=(1/2*solver.nuec-solver.gammaj);

solver.Ek2=zeros(size(solver.nuej));

if 1
	allk=[min(k) max(k)];
	subplot(211), plot(k, solver.Femj, k, solver.Febj, k, solver.Fej)
	set(gca,'xlim',allk); title('F_e(k)');
	subplot(212), plot(k, solver.nuej,k,solver.nue)
	set(gca,'xlim',allk); title('\nu_e(k)');
	%pause
end

% Initialise with numerically estimated from electron velocity distribution
solver.nue=solver.nuej;

% Boundary conditions as in 
% Sanbonmatsu et al., Quantitative comparison of reduced-description
% particle-in-cell and quasilinear-Zakharov models for parametrically
% excited Langmuir turbulence, Phys. Plasma, 7, 2824-2841, 2000
solver.ii=find(solver.k~=0.0 & ...
  abs(solver.vj)>1.3*solver.ve & abs(solver.vj)<333*solver.ve);
solver.ip=find(solver.vj>1.3*solver.ve & solver.vj<333*solver.ve);
solver.im=find(solver.vj<-1.3*solver.ve & solver.vj>-333*solver.ve);

function dy=derivate(dx,y)
% function dy=derivate(dx,y)
dy = zeros(size(y(:)));

dy(1) = (y(2)-y(1))/dx;
dy(2:end-1) = (y(3:end)-y(1:end-2))/(2.0*dx);
dy(end) = (y(end)-y(end-1))/dx;

