function solver=collisions(solver)
% function solver=collisions(solver)

%
% $Id: collisions.m,v 1.10 2011/03/26 09:20:41 patrick Exp $
%
% Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
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

ik0 = find(solver.k==0);
ik = find(solver.k~=0);

% Ion Landau damping of ion acoustic waves 
% from Hanssen et al., Numerical test of the weak turbulence approximation
% to ionospheric Langmuir turbulence, JGR, 97, 12073-12091, 1992 
if ~isfield(solver,'nui0'),
	solver.nui0=sqrt(pi/8)*(sqrt(si.me/solver.mi)+ ...
		(solver.Te/solver.Ti)^1.5*exp(-solver.Te/(2*solver.Ti)-3/2));
end

% Ion damping ( collisions + Landau damping)
% from Hanssen et al., Numerical test of the weak turbulence approximation
% to ionospheric Langmuir turbulence, JGR, 97, 12073-12091, 1992 
solver.nui=solver.nuic/solver.tau+solver.nui0*abs(solver.k);

%solver.nui=zeros(size(solver.nui));

if isfield(solver,'lhanssen') & solver.lhanssen
% Electron background Landau damping
% from Hanssen et al., Numerical test of the weak turbulence approximation
% to ionospheric Langmuir turbulence, JGR, 97, 12073-12091, 1992 
	solver.gammal=zeros(size(solver.k));
	solver.gammal(ik)=sqrt(pi/8)*9/4*(solver.mi/(solver.eta*si.me))^2.5./ ...
		abs(solver.k(ik)).^3.*...
		exp(-9/8*solver.mi/(solver.eta*si.me)./(solver.k(ik).^2)-3/2);
else
% Electron background Landau damping calculated from the velocity distribution
% function of the electrons
	k=solver.k/solver.chi;
	wpe=solver.wpe;
	ne=solver.ne;
	ve=solver.ve;
	ve2=solver.ve2;
	tau=solver.tau;
	ke=1/solver.lambdae;
	Fe=zeros(size(k));
	Fe(ik)=ne/sqrt(2*pi)/ve*exp(-(ke./k(ik)).^2/2);
	dFe=zeros(size(k));
	dFe(ik)=Fe(ik).*(-wpe./k(ik)/ve2);
	solver.gammal=zeros(size(k));
	%solver.gammal(ik)=-pi/2/ne*wpe^3./k(ik).^2.*sign(k(ik)).*dFe(ik)*0.0992*tau;
	solver.gammal(ik)=-pi/2/ne*wpe^3./k(ik).^2.*sign(k(ik)).*dFe(ik)*tau;
end

% Beam growth rate
if isfield(solver,'k0'),
	solver.gammab=zeros(size(solver.k));
if 0
	u=solver.k/solver.k0;
	solver.gammab(ik)=solver.gamma0*(u(ik)-1)./((solver.dvb_vb*u(ik)).^3).* ...
		exp(-((u(ik)-1)./(solver.dvb_vb*u(ik))).^2);
	subplot(211), plot(solver.k,solver.gammab);
end
	u=zeros(size(solver.k));
	u(ik)=solver.k0./solver.k(ik);
	solver.gammab(ik)= solver.gamma0 * ...
		(u(ik)/solver.dvb_vb).^2.*((1-u(ik))./solver.dvb_vb).* ...
		exp(-((1-u(ik))./solver.dvb_vb).^2);
	solver.gammab(ik0)=0;
if 0
	subplot(212), plot(solver.k,solver.gammab);
	keyboard
end
else
	solver.gammab=zeros(size(solver.k));
end

% Electron damping (Landau + collision + beam)
solver.nue=(1/2*solver.nuec+solver.gammal-solver.gammab);

solver.nueL=(1/2*solver.nuec+solver.gammal);


