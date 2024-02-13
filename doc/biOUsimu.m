function biOUsimu
% function biOUsimu

%
% $Id: biOUsimu.m,v 1.8 2010/06/17 17:26:22 patrick Exp $
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

seed = 10;
nexp = 300;

nu = 1;
k  = 3;

s  = 2;

D  = sqrt(nu*s);
X0 = sqrt(s)*exp(i*2*pi*.2);;

if 0
h  = 1e-3;
n   = 6000;
cmp(X0, h, n, nu, k, s, nexp, seed)
figure
h  = 1e-4;
n   = 6000;
cmp(X0, h, n, nu, k, s, nexp, seed)
end

h  = 1e-3;
n   = 5000;
cmp(X0, h, n, nu, k, s, nexp, seed)

return
orient landscape; set(gcf,'PaperOrientation','portrait');
exportfig(gcf,'biOUsimu1','color','cmyk','boundscode','mcode');

figure
h  = 1.6;
n   = fix(1000/h);
cmp(X0, h, n, nu, k, s, nexp, seed)
orient landscape; set(gcf,'PaperOrientation','portrait');
exportfig(gcf,'biOUsimu2','color','cmyk','boundscode','mcode');

figure
h  = 1.99;
n   = fix(1000/h);
cmp(X0, h, n, nu, k, s, nexp, seed)
orient landscape; set(gcf,'PaperOrientation','portrait');
exportfig(gcf,'biOUsimu3','color','cmyk','boundscode','mcode');


function cmp(X0, h, n, nu, k, s, nexp, seed)

if 0,
% normalise by standard deviation
scale = sqrt(s*nu);
else
scale = 1;
end

RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed));
Se = biOU_exact(X0, h, n, nu, k, s, nexp);
Se = Se./scale;
Sem = mean(Se).';
Ses = std(Se).';
Secov = mean(real(Se).*imag(Se))';

RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed));
Sa = biOU_euler(X0, h, n, nu, k, s, nexp);
Sa = Sa./scale;
Sam = mean(Sa)';
Sas = std(Sa)';

t = h*[0:n]';
Sm = biOUmean(t,X0,nu,k); 
[Svar,Scov] = biOUvar(t,X0,nu,k,s);
Scov = Scov + real(Sm).*imag(Sm);
Ss = biOUsdev(t,X0,nu,k,s);
Sm = Sm./scale;
Ss = Ss./scale;

subplot(311), plot(t,[real(Sm),imag(Sm),abs(Sm)]); legend('Re','Im','abs')
subplot(312), plot(t,[real(Svar),imag(Svar),abs(Scov)]); legend('Re','Im','cov')
subplot(313), plot(t,[real(Ss),imag(Ss),abs(Ss)]); legend('Re','Im','abs')
pause

%plot(t,[Se',Sa'],t,[Sm,Sm+Ss,Sm-Ss]);

specs = sprintf('%d paths @ h=%g, \\nu=%g k=%g \\sigma=%g', nexp, h, nu, k, s);
tlim = [min(t), max(t)];


subplot(211),
plot(...t,mean(abs(Sa)),'r',... %t,imag(Sa(1,:)),'r',t,mean(abs(Sa)),...
     t,real([Sam,Sam+Sas,Sam-Sas]),'b',...%t,imag([Sam,Sam+Sas,Sam-Sas]),'b',...
     t,real([Sm,Sm+Ss,Sm-Ss]),'k')%,t,imag([Sm,Sm+Ss,Sm-Ss]),'k');
set(gca,'xlim',tlim);
ylim = get(gca,'ylim');
title(sprintf('Euler-Maruyama scheme (%s)',specs))
if scale==1,
ylabel('X(t)');
else
ylabel('X(t)/(c\tau/2)^{1/2}');
end

subplot(212),
plot(...t,mean(abs(Se)),'r',...%t,imag(Se(1,:)),'r',t,mean(abs(Se)),'r',...
     t,real([Sem,Sem+Ses,Sem-Ses]),'b',...%t,imag([Sem,Sem+Ses,Sem-Ses]),'b',...
     t,real([Sm,Sm+Ss,Sm-Ss]),'k')%,t,imag([Sm,Sm+Ss,Sm-Ss]),'k');
set(gca,'xlim',tlim);
set(gca,'ylim',ylim);
title(sprintf('Exact scheme (%s)',specs))
if scale==1,
ylabel('X(t)');
else
ylabel('X(t)/(c\tau/2)^{1/2}');
end


function X = biOU_exact(X0, h, n, nu, k, s, nexp)

if 0
Nt = 1/2*sqrt(s^2/nu*(1-exp(-2*nu*h)))*(randn(nexp,n)+i*randn(nexp,n));
else
[v,c] = biOUvar(h,X0,nu,k,s);
if 1
N1 = randn(nexp,n);
N2 = randn(nexp,n);
Nt = sqrt(real(v))*N1+i*...
     (sqrt(imag(v)-c/real(v))*N2+sqrt(c/real(v))*N1);

else
Nt = sqrt(real(v/2)).*randn(nexp,n)+sqrt(imag(v/2)).*i*randn(nexp,n);
end
end
X = zeros(nexp,n+1);

a = exp(-(i*k^2+nu)*h);

%S = [X0 filter(1,[1, -a], s*Nt, a*S0)];

% trick to use filter to perform the recursion
% S(0) = S0
% S(n+1) = S(n)*a +  b(n)
for k=1:nexp,
  X(k,:) = filter(1,[1, -a], [X0, Nt(k,:)]);
end


function X = biOU_euler(X0, h, n, nu, k, s, nexp)

Nt = sqrt(s^2/(2*nu)*h)*(randn(nexp,n)+i*randn(nexp,n));
X = zeros(nexp,n+1);
a = exp(-(i*k^2+nu)*h);
for k=1:nexp,
  X(k,:) = filter(1,[1, -a], [X0, Nt(k,:)]);
end

