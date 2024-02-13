function solver=test_cerenkov_diffusion
% function solver=test_cerenkov_diffusion

%
% $Id: test_cerenkov_diffusion.m,v 1.15 2011/03/26 09:20:42 patrick Exp $
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

% spatial period
solver.L=15*pi;
% spatial grid size
solver.N=4096;
% time increment
solver.h=1e-3;

% initial E-field power
solver.E0=0.0; 
% initial E-field fluctuation
solver.Erms=1e-4; 

solver.edrive=0.0; % electric pump field
solver.domega=0.0; % frequency mismatch

% damping
solver.nuic=0.01; % ion collisions [s^-1]
solver.nuec=2.4;  % electron collisions [s^-1]
solver.nui0=0.3;  % dimensionless ion Landau damping 

% temperatures
solver.Te=3000.0;
solver.Ti=1200.0;

% beam parameter
solver.gamma0=2.0*5;
%solver.k0=-22.0;
solver.k0=-10.0;
solver.dvb_vb=0.3;

% plasma parameters
solver.ne=7.e10;
solver.mui=16; % Molecular oxygen number

% si k0 est libre decomenter la ligne suivante 
if ~isfield(solver,'k0'), solver.vb=6.6235e6; end

% Initialise background parameters
solver = maxwell_param_init(solver);
% Initialise beam parameters
solver = beam_param_init(solver);
% Collision frequencies
solver = collisions(solver);
% Initialise thermal fluctuations
solver = thermal_fluctuations_init(solver);
% Initialise diffusion
solver = diffusion_init(solver);

solver.display=1;
solver.mdisplay=100;

solver.save=1;
solver.msave=100;

solver.maxiter = 20480;
solver.maxtime = solver.maxiter*solver.h;

solver.av_start=1;
solver.av_end=solver.maxiter;

% Wave vectors to consider to calculate IS spectra
% Radar parameters
si=SIconstants;
f=224e6; 
solver.is_kmin=2*2*pi*f/si.c; % [m-1]
f=931e6;
solver.is_kmax=2*2*pi*f/si.c; % [m-1]
solver.is_nbk=12;

solver.is_start=1;
solver.is_end=solver.maxiter;
solver.is_stride=1;
%solver.is_nfft=1024;

solver.ctime=0;
solver.times=[];
solver.Fejs=[];
solver.dFejs=[];
for i=1:solver.maxiter,

  solver.citer=i;
	solver.ctime = solver.ctime + solver.h;

	solver = diffusion_update(solver);

	if mod(i,100)==0,
		s=sprintf(' Iter %6.d/%6.d: integrated to %.2f/%.2f (%.2f%%)',i, ...
			solver.maxiter,solver.ctime,solver.maxtime, ...
			100*solver.ctime/solver.maxtime);
		fprintf(1,'%s%s', s, repmat(8,1,length(s)));
	end

	if solver.save & mod(i,solver.msave)==0,
		solver.times=[solver.times solver.ctime];
		solver.Fejs=[solver.Fejs, solver.Fej];
		solver.dFejs=[solver.dFejs, solver.Femj-solver.Fej];
	end

end;
if solver.maxiter,
	fprintf(1,' Iter %6.d/%6.d: integrated to %.2f/%.2f (%.2f%%)\n',i, ...
		solver.maxiter,solver.ctime,solver.maxtime, ...
		100*solver.ctime/solver.maxtime);
end

function solver=diffusion_update(solver)
% function solver=diffusion_update(solver)

si=SIconstants;

ip=solver.ip;
im=solver.im;

k=solver.k/solver.chi;
L=solver.L*solver.chi;
ne=solver.ne;
wpe=solver.wpe;
tau=solver.tau;
chi=solver.chi;
dt=solver.h*solver.tau;
v=solver.vj;
ve2=solver.ve2;

% Calculate D(v,t) from Eq. 10 of
% Sanbonmatsu et al., Quantitative comparison of reduced-description
% particle-in-cell and quasilinear-Zakharov models for parametrically
% excited Langmuir turbulence, Phys. Plasma, 7, 2824-2841, 2000
%C=pi*wpe*eps0/ne/me/L;
C=4*pi^2*si.e^2/si.me^2/L;

dD=C*abs(k/wpe).*abs(solver.delta_E).^2*solver.epsilon^2;

dD1=zeros(size(dD));
nue=solver.gammal/tau;
ii=find(nue~=0.0 & (k.*v-wpe)~=0.0);
for i=1:solver.N,
	term=abs(solver.delta_E(ii)).^2*solver.epsilon^2.*nue(ii)./ ...
		((k(i).*v(ii)-wpe).^2+nue(ii).^2);
	dD1(i)=e^2/me^2*sum(term);
end
subplot(211), plot(k,dD)
subplot(212), plot(k,dD1)
keyboard
pause

if rem(solver.citer,10)==0,
	if isfield(solver,'vb'),
		xlim=solver.wpe/solver.vb*[1-sign(solver.vb)*0.8 1+sign(solver.vb)*0.8];
	else
		xlim=solver.wpe/solver.ve*[0.5 2];
	end
	allk=[min(k) max(k)];
	subplot(511), plot(k(im), dD(im), k(ip), dD(ip)), set(gca,'xlim',allk);
	subplot(512), plot(k, solver.Fej), set(gca,'xlim',xlim);
	Fej=solver.Fej;
end

dk=solver.dk;
dk_dv=solver.dk_dv;

solver.Fej(ip) = ftcs_lax(dk, dt, solver.Fej(ip), v(ip).*dD(ip)/ve2, dk_dv(ip), solver);
solver.Fej(im) = ftcs_lax(dk, dt, solver.Fej(im), v(im).*dD(im)/ve2, dk_dv(im), solver);

solver.dFej(ip)=derivate(dk,solver.Fej(ip)).*dk_dv(ip);
solver.dFej(im)=derivate(dk,solver.Fej(ip)).*dk_dv(im);

%solver.gammakj=zeros(size(solver.Ek));
solver.gammakj=solver.gammaj;
solver.gammakj(ip)=-pi/2/ne*wpe^3./k(ip).^2.*sign(k(ip)).*solver.dFej(ip)*tau;
solver.gammakj(im)=-pi/2/ne*wpe^3./k(im).^2.*sign(k(im)).*solver.dFej(im)*tau;
solver.nue=(1/2*solver.nuec+solver.gammakj);

if rem(solver.citer,10)==0,
	subplot(513), plot(k, solver.Fej-Fej), set(gca,'xlim',xlim);
	subplot(514), plot(k, solver.Fej), set(gca,'xlim',xlim);
	subplot(515), plot(k, solver.nue), set(gca,'xlim',xlim);
	drawnow
	%pause
end
%subplot(212), plot(k, solver.nue), drawnow, keyboard

function u=ftcs_lax(dx,dt,u,D,dk_dv,solver)
% FTCS corrected by Lax method 

% Courant condition for stability
c=abs(dk_dv).*D*dt/dx^2;
if max(c)>1, fprintf('Error: max(c)>1\n'); keyboard; end

u(2:end-1) = 1/2*(u(1:end-2)+u(3:end)) + dk_dv(2:end-1)*dt/(2*dx).* ...
	(D(3:end).*u(3:end)-D(1:end-2).*u(1:end-2));  


function dy=derivate(dx,y)
% function dy=derivate(dx,y)
dy = zeros(size(y(:)));

dy(1) = (y(2)-y(1))/dx;
dy(2:end-1) = (y(3:end)-y(1:end-2))/(2*dx);
dy(end) = (y(end)-y(end-1))/dx;

