function s=zak0(paramfile)
% function s=zak0(paramfile)

%
% $Id: zak0.m,v 1.15 2011/03/26 09:20:42 patrick Exp $
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
	 ftdn0 = fft(dn0);
elseif s.choice == 2,
end

Ek2s = zeros(s.n,s.ns);
nk2s = zeros(s.n,s.ns);
nEm = zeros(s.ns,1);

tic
for count=0:s.nq-1,
  Ek2 = fftshift(abs(fft(e0)).^2);
	nk2 = fftshift(abs(fft(n0)).^2);
  if rem(count,100)==0,
	  fprintf(1,'iter %8d/%8d%s', count, s.nq, char(8*ones(1,23)))
		if 0
		fprintf(1,'%10e %10e\t', Ek2([s.klim(1),s.klim(1)+1]));
		fprintf(1,'%10e %10e\n', Ek2([s.klim(2)-1,s.klim(2)]));
		fprintf(1,'%10e %10e\t', nk2([s.klim(1),s.klim(1)+1]));
		fprintf(1,'%10e %10e\n', nk2([s.klim(2)-1,s.klim(2)]));
		end
		if 1
		subplot(211), plot(Ek2), set(gca,'xlim',s.klim);
		subplot(212), plot(nk2), set(gca,'xlim',s.klim);
		drawnow
		end
	end

  jk = fix(count/fix(s.nq/s.ns))+1;
	Ek2s(:,jk) = Ek2s(:,jk) + Ek2;
	nk2s(:,jk) = nk2s(:,jk) + nk2;
	nEm(jk) = nEm(jk) + mean(n0.*e0);
	if 0 & rem(count,fix(s.nq/s.ns))==fix(s.nq/s.ns)-1
		subplot(211), imagesc(Ek2s'), axis xy,set(gca,'xlim',s.klim); 
		subplot(212), imagesc(nk2s'), axis xy,set(gca,'xlim',s.klim);
		drawnow
	end

if 0
  % ponderomotive mismatch
	% long term correction due to density modification
	% introduced by Harvey Rose
  s.w = sum(abs(fte0(2:end))^2); 
else
  s.w = 0.0;
end

	fte02 = fft(abs(s.edrive+e0).^2);
	ftn0 = fft(n0);

	ftn1 = ftn1sub(ftn0,ftdn0,fte02,s);
	n1 = ifft(ftn1);

	e1 = e1sub(n1,n0,e0,s);

	ftdn0 = ftdn0sub(ftdn0,fte02,ftn1,ftn0,e1,s);

	e0 = e1;
	n0 = n1;
end
trun=fix(toc);
fprintf(1,'\nCPU time %2d:%2d [min:sec]\n', fix(trun/60), mod(trun,60));

s.Ek2s = Ek2s / fix(s.nq/s.ns);
s.nk2s = nk2s / fix(s.nq/s.ns);
s.nEm = nEm / fix(s.nq/s.ns);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ftn1=ftn1sub(ftn0,ftdn0,fte02,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = [1:s.mmax]'; 

k2 = s.pi2ndx2*m.^2;
nyik = s.nyi(m);
nyik2 = nyik.^2;
e = exp(-nyik*s.h);
sq = sqrt(k2-nyik2);
frac  = zeros(size(m));
frac(find(sq == 0.0)) = 1;
ii = find(sq ~= 0.0);
frac(ii) = sin(sq(ii)*s.h)./sq(ii);
co = cos(sq*s.h);

ftn1 = zeros(size(ftn0));

ftn1(m+1)     = e.*frac.*(nyik.*ftn0(m+1)+ftdn0(m+1) ...
                -s.h2*k2.*fte02(m+1))+e.*co.*ftn0(m+1);
ftn1(s.n+1-m) = e.*frac.*(nyik.*ftn0(s.n+1-m)+ftdn0(s.n+1-m) ...
                -s.h2*k2.*fte02(s.n+1-m))+e.*co.*ftn0(s.n+1-m);

% k=0 mode
ftn1(1) = ftn0(1)+ftdn0(1);

if m(end) == s.n2-1, % k=-(n/2)dk mode
  m = s.n2+1;
  k2 = s.pi2ndx2*s.n2^2;
  nyik = s.nyi(m);
  nyik2 = nyik^2;
  e = exp(-nyik*s.h);
  sq = sqrt(k2-nyik2);
  if sq == 0.0, 
    frac = 1;
  else
    frac = sin(sq*s.h)./sq;
  end
  co = cos(sq*s.h);
  ftn1(m) = e*frac*(nyik*ftn0(m)+ftdn0(m)-s.h2*k2*fte02(m)) ...
            +e*co*ftn0(m);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e1=e1sub(n1,n0,e0,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imh = s.im*s.h;
c = imh/2;

ne0 = n0.*(s.edrive+e0);
ftne0 = fft(ne0);
fte0 = fft(e0);

ftf = zeros(size(n1));
ftg = zeros(size(n1));

m = [1:s.mmax]'; % 2*mmax+1 modes calculated
e = exp(-(imh*(s.pi2ndx2*m.^2-s.domega-s.w) + s.nye(m)*s.h));
ftf(m+1)     = e.*fte0(m+1);
ftf(s.n+1-m) = e.*fte0(s.n+1-m);
ftg(m+1)     = e.*ftne0(m+1);
ftg(s.n+1-m) = e.*ftne0(s.n+1-m);

% k=0 mode
e = exp(imh*(s.domega+s.w)-s.nye(1)*s.h);
n0e0m = mean(n0.*e0); % source term <nE>
ftf(1) = e*fte0(1);
ftg(1) = e*(ftne0(1)-n0e0m);

if m(end) == s.n2-1, % k=-(n/2)dk mode
  m = s.n2+1;
  k2 = s.pi2ndx2*s.n2^2;
  e = exp(imh*(k2-s.domega-s.w) + s.nye(m)*s.h);
  ftf(m) = e*fte0(m);
  ftg(m) = e*ftne0(m);
end

f = ifft(ftf);
g = ifft(ftg);

n1e0m = mean(n1.*e0); % source term <nE>

e1 = (f - c*(g + s.edrive*n1 - n1e0m))./(1+c*n1);

% zero padding of aliased modes when k=0 is first element
fte1 = fft(e1);
fte1(s.zm) = 0;
e1 = ifft(fte1);

if 0
Ek2 = fftshift(abs(fft(e1)).^2);
fprintf(1,'%10e %10e\t', Ek2([s.klim(1),s.klim(1)+1]));
fprintf(1,'%10e %10e\n', Ek2([s.klim(2)-1,s.klim(2)]));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ftdn0 = ftdn0sub(ftdn0,fte02,ftn1,ftn0,e1,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e12 = abs(s.edrive+e1).^2;

fte12 = fft(e12);

m = [1:s.mmax]';

e = exp(-2*s.h*s.nyi(m));
c = s.h2*s.pi2ndx2*m.^2;
ftdn0(m+1)     = e.*(ftdn0(m+1)-c.*(ftn0(m+1)+fte02(m+1))) ...
                 -c.*(ftn1(m+1)+fte12(m+1));
ftdn0(s.n+1-m) = e.*(ftdn0(s.n+1-m)-c.*(ftn0(s.n+1-m)+fte02(s.n+1-m))) ...
	               -c.*(ftn1(s.n+1-m)+fte12(s.n+1-m));

% k=0 mode
ftdn0(1) = exp(-2*s.h*s.nyi(1))*ftdn0(1); 

if m(end) == s.n2-1, % k=(-n/2)dk  mode
  m = s.n2+1; 
  e = exp(-2*s.h*s.nyi(m));
  c = s.h2*s.pi2ndx2*s.n2^2;
  ftdn0(m+1) = e*(ftdn0(m+1)-c*(ftn0(m+1)+fte02(m+1)))-c*(ftn1(m+1)+fte12(m+1));
end


