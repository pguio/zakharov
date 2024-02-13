function s=zak2(paramfile)
% function s=zak2(paramfile)

%
% $Id: zak2.m,v 1.10 2011/03/26 09:20:42 patrick Exp $
%
% Copyright (c) 2001-2011 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

s = feval(paramfile);
s = nyisub(s);
s = nyesub(s);

if s.choice == 1,
   [n0,dn0,e0] = initial(s);
	 ftn0  = fft(n0);
	 ftdn0 = fft(dn0);
	 fte0  = fft(e0);
elseif s.choice == 2,
end

Ek2s = zeros(s.n,s.ns);
nk2s = zeros(s.n,s.ns);
buff = zeros(s.n,1);
nEm = zeros(s.ns,1);

options=odeset('RelTol',1e-3,'Stats','on');

tic
for count=0:s.odesteps:s.nq-1,
  Ek2 = fftshift(abs(fte0).^2);
	nk2 = fftshift(abs(ftn0).^2);
  if rem(count,10)==0,
	  fprintf(1,'iter %8d/%8d\n', count, s.nq)
		if 1
		subplot(211), plot(Ek2), set(gca,'xlim',s.klim);
		subplot(212), plot(nk2), set(gca,'xlim',s.klim);
		drawnow
		end
	end

	if 0 & rem(count,fix(s.nq/s.ns))==fix(s.nq/s.ns)-1
		subplot(211), imagesc(Ek2s'),axis xy,set(gca,'xlim',s.klim); 
		subplot(212), imagesc(nk2s'),axis xy,set(gca,'xlim',s.klim);
		drawnow
	end

	Xk0 = [fte0(s.cm), ftn0(s.cm), ftdn0(s.cm)];
	tspan = [0:s.odesteps]*s.h;

% fprintf(1,'%f %f\n',s.t0+[tspan(1),tspan(end)]);
	[t,Xk1]=ode45(@zakderiv,tspan,Xk0,options,s);
	s.t0 = s.t0 + tspan(end);

	fte0(s.cm) = Xk1(end,1:s.ncm);
	ftn0(s.cm) = Xk1(end,s.ncm+1:2*s.ncm);
	ftdn0(s.cm) = Xk1(end,2*s.ncm+1:3*s.ncm);

  jk = fix(count/fix(s.nq/s.ns))+1;
  buff(s.cm) = sum(abs(Xk1(1:end-1,1:s.ncm)).^2,1);
  Ek2s(:,jk) = Ek2s(:,jk)+fftshift(buff,1);
  buff(s.cm) = sum(abs(Xk1(1:end-1,s.ncm+1:2*s.ncm)).^2,1);
  nk2s(:,jk) = nk2s(:,jk)+fftshift(buff,1);

	tE = zeros(s.n,s.odesteps);
	tE(s.cm,:) = Xk1(1:end-1,1:s.ncm).';
	tn = zeros(s.n,s.odesteps);
	tn(s.cm,:) = Xk1(1:end-1,s.ncm+1:2*s.ncm).';
	nEm(jk) = nEm(jk) + sum(sum(ifft(tn).*ifft(tE)))/s.n;

end
trun=fix(toc);
fprintf(1,'\nCPU time %2d:%2d [min:sec]\n', fix(trun/60), mod(trun,60));

s.Ek2s = Ek2s / fix(s.nq/s.ns);
s.nk2s = nk2s / fix(s.nq/s.ns);
s.nEm = nEm / fix(s.nq/s.ns);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dXk = zakderiv(t,Xk,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ek = zeros(s.n,1);
nk = zeros(s.n,1);
vk = zeros(s.n,1);

Ek(s.cm) = Xk([1:s.ncm]');
nk(s.cm) = Xk([s.ncm+1:2*s.ncm]');
vk(s.cm) = Xk([2*s.ncm+1:3*s.ncm]');

if 0
  % ponderomotive mismatch
  % long term correction due to density modification
  % introduced by Harvey Rose
  s.w = sum(abs(fte0(2:end))^2);
else
  s.w = 0.0;
end

Ej = ifft(Ek);
nj = ifft(nk);

njEjm = mean(nj.*Ej); % source term <nE>

Ejtot = Ej+s.edrive*exp(-i*(s.domega+s.w)*(s.t0+t));

nEk = fft(nj.*Ejtot);

% zero padding of aliased modes when k=0 is first element
nEk(s.zm) = 0;

j = s.cm;

dEk = -(s.nye(j) + s.ik2(j)).*Ek(j) -i*nEk(j);
dEk(1) = dEk(1) + i*njEjm;

dnk = vk(j);

E2k = fft(abs(Ejtot).^2);
dvk = -2*s.nyi(j).*vk(j) - s.k2(j).*(nk(j) + E2k(j));

dXk = [dEk; dnk; dvk];

