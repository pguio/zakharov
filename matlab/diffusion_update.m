function solver=diffusion_update(solver)
% function solver=diffusion_update(solver)

%
% $Id: diffusion_update.m,v 1.31 2011/03/26 09:20:41 patrick Exp $
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

% Boundary conditions as in 
% Sanbonmatsu et al., Quantitative comparison of reduced-description
% particle-in-cell and quasilinear-Zakharov models for parametrically
% excited Langmuir turbulence, Phys. Plasma, 7, 2824-2841, 2000

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
C=4*pi^2*si.e^2/si.me^2/L;
%D=C*abs(k/wpe).*abs(solver.Ek/solver.N*solver.epsilon).^2;
D=C*abs(k/wpe).*abs(solver.Ek*solver.epsilon).^2;

dD=C*abs(k/wpe).*abs(solver.delta_Ek/solver.N*solver.epsilon).^2;
vdD_ve2=v.*dD/ve2;

dk=solver.dk;
dk_dv=solver.dk_dv;

solver.Fej(im) = cn2(dk,dt,solver.Fej(im),D(im),dk_dv(im),solver);
solver.Fej(ip) = cn2(dk,dt,solver.Fej(ip),D(ip),dk_dv(ip),solver);
%solver.Fej(im) = ftcs_lax(dk,dt,solver.Fej(im),vdD_ve2(im),dk_dv(im),solver);
%solver.Fej(ip) = ftcs_lax(dk,dt,solver.Fej(ip),vdD_ve2(ip),dk_dv(ip),solver);

solver.dFej(im)=derivate(dk,solver.Fej(im)).*dk_dv(im);
solver.dFej(ip)=derivate(dk,solver.Fej(ip)).*dk_dv(ip);

solver.gammakj=solver.gammaj;
solver.gammakj(im)=pi/2/ne*wpe^3./k(im).^2.*sign(k(im)).*solver.dFej(im)*tau;
solver.gammakj(ip)=pi/2/ne*wpe^3./k(ip).^2.*sign(k(ip)).*solver.dFej(ip)*tau;
solver.nue=(1/2*solver.nuec-solver.gammakj);


function u=ftcs_lax(dx,dt,u,D,dk_dv,solver)
% FTCS scheme corrected by Lax method 

% Courant condition for stability
c=abs(dk_dv).*D*dt/dx^2;
if max(c)>1, 
	fprintf('Error: max(c)>1\n'); 
	keyboard; 
end

a=dt/(2*dx^2);

u(2:end-1) = 1/2*(u(1:end-2)+u(3:end)) + dk_dv(2:end-1)*a.* ...
	(D(3:end).*u(3:end)-D(1:end-2).*u(1:end-2));  


function u=ftcs2(dx,dt,u,D,dk_dv,solver)
% FTCS scheme fully explicit
% with Dirichlet boundary condition

% Courant like condition for stability
c=2*D.*dk_dv*dt/dx^2;
if max(c)>1, 
	fprintf('Error: max(c)>1\n'); 
	keyboard; 
end

a=dt/(2*dx^2);

Di = (D(1:end-1)+D(2:end))/2;
dk_dvi = (dk_dv(1:end-1)+dk_dv(2:end))/2;
Di = Di.*dk_dvi;
u(2:end-1) = u(2:end-1)+dk_dv(2:end-1)*a.* ...
	(Di(2:end).*(u(3:end)-u(2:end-1))-Di(1:end-1).*(u(2:end-1)-u(1:end-2)));

function u=cn2(dx,dt,u,D,dk_dv,solver)
% Crank-Nicholson scheme (semi-implicit method)
% with Dirichlet boundary condition

n=length(D);

Di = (D(1:end-1)+D(2:end))/2;
dk_dvi = (dk_dv(1:end-1)+dk_dv(2:end))/2;
Di = Di.*dk_dvi;
a=dt/(2*dx^2);
A = spdiags( ...
	[[-dk_dv(2:end-1)*a.*Di(1:end-1);0;0], ...
	[1;1+dk_dv(2:end-1)*a.*(Di(1:end-1)+Di(2:end));1], ...
	[0;0;-dk_dv(2:end-1)*a.*Di(2:end)]], [-1:1], n, n);
b = [u(1); ...
	(dk_dv(2:end-1)*a.*Di(1:end-1).*u(1:end-2)+...
	(1-dk_dv(2:end-1)*a.*(Di(1:end-1)+Di(2:end))).*u(2:end-1)+...
	dk_dv(2:end-1)*a.*Di(2:end).*u(3:end)); u(end)];

[c,v]=condest(A);
if norm(u,1)/c<eps,
	keyboard
	[L,U]=luinc(A,1e-6);
%	u = bicgstab(A,b,1e-15,10,L,U);
%	u = bicg(A,b,1e-15,10,L,U);
%	u = cgs(A,b,1e-15,10,L,U);
	u = gmres(A,b,1e-15,10,L,U);
%	u = qmr(A,b,1e-15,10,L,U);
%	R = cholinc(A,1e-3);
%	u = pcg(A,b,1e-12,10,R',R);
else
	u = A \ b;
end

function dy=derivate(dx,y)
% function dy=derivate(dx,y)

dy = zeros(size(y(:)));

dy(1) = (y(2)-y(1))/dx;
dy(2:end-1) = (y(3:end)-y(1:end-2))./(2.0*dx);
dy(end) = (y(end)-y(end-1))/dx;

