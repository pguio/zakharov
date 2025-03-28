function OUsimu
% function OUsimu

%
% $Id: OUsimu.m,v 1.6 2025/03/24 15:53:53 patrick Exp $
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

c   = 1;
tau = 1;

mu  = 0;
%S0  = sqrt(c*tau/2);
S0 = 0;

if 0
dT  = 1e-3;
n   = 6000;
cmp(S0, dT, n, mu, c, tau, nexp, seed)
figure
dT  = 1e-4;
n   = 6000;
cmp(S0, dT, n, mu, c, tau, nexp, seed)
figure
c   = 1e6;
tau = 1e-3;
cmp(S0, dT, n, mu, c, tau, nexp, seed)
end

dT  = 1e-3;
n   = 5000;
cmp(S0, dT, n, mu, c, tau, nexp, seed)
orient landscape; set(gcf,'PaperOrientation','portrait');
exportfig(gcf,'OUsimu1','color','cmyk','boundscode','mcode');

figure
dT  = 1.6;
n   = fix(1000/dT);
cmp(S0, dT, n, mu, c, tau, nexp, seed)
orient landscape; set(gcf,'PaperOrientation','portrait');
exportfig(gcf,'OUsimu2','color','cmyk','boundscode','mcode');

figure
dT  = 1.99;
n   = fix(1000/dT);
cmp(S0, dT, n, mu, c, tau, nexp, seed)
orient landscape; set(gcf,'PaperOrientation','portrait');
exportfig(gcf,'OUsimu3','color','cmyk','boundscode','mcode');


function cmp(S0, dT, n, mu, c, tau, nexp, seed)

if 1,
% normalise by standard deviation
scale = sqrt(c*tau/2);
else
scale = 1;
end

RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));
Se = OU_exact(S0, dT, n, mu, c, tau, nexp);
Se = Se./scale;
Sem = mean(Se)';
Ses = std(Se)';

RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));
Sa = OU_euler(S0, dT, n, mu, c, tau, nexp);
Sa = Sa./scale;
Sam = mean(Sa)';
Sas = std(Sa)';

t = dT*[0:n]';
Sm = mu+(S0-mu)*exp(-t/tau);
Ss = sqrt(c*tau/2*(1-exp(-2*t/tau)));
Sm = Sm./scale;
Ss = Ss./scale;

%plot(t,[Se',Sa'],t,[Sm,Sm+Ss,Sm-Ss]);

specs = sprintf('%d paths @ h=%g, \\tau=%g c=%g', nexp, dT, tau, c);
tlim = [min(t), max(t)];


subplot(211),
plot(t,Sa(1,:),'r',t,[Sam,Sam+Sas,Sam-Sas],'b',t,[Sm,Sm+Ss,Sm-Ss],'k');
set(gca,'xlim',tlim);
ylim = get(gca,'ylim');
title(sprintf('Euler-Maruyama scheme (%s)',specs))
if scale==1,
ylabel('X(t)');
else
ylabel('X(t)/(c\tau/2)^{1/2}');
end

subplot(212),
plot(t,Se(1,:),'r',t,[Sem,Sem+Ses,Sem-Ses],'b',t,[Sm,Sm+Ss,Sm-Ss],'k');
set(gca,'xlim',tlim);
set(gca,'ylim',ylim);
title(sprintf('Exact scheme (%s)',specs))
if scale==1,
ylabel('X(t)');
else
ylabel('X(t)/(c\tau/2)^{1/2}');
end


function S = OU_exact(S0, dT, n, mu, c, tau, nexp)

a = exp(-dT/tau);
m = mu*(1-a);
s = sqrt(c*tau/2*(1-a*a));

fprintf(1,'a=%.6f m=%.6f s=%.6f\n',a,m,s);

Nt = randn(nexp,n);
S = zeros(nexp,n+1);
%S = [S0 filter(1,[1, -a], m+s*Nt, a*S0)];

% trick to use filter to perform the recursion
% S(0) = S0
% S(n+1) = S(n)*a +  b(n)
for i=1:nexp,
  S(i,:) = filter(1,[1, -a], [S0, m+s*Nt(i,:)]);
end


function S = OU_euler(S0, dT, n, mu, c, tau, nexp)

a = 1 - dT/tau;
m = mu*dT/tau;
s = sqrt(c*dT);

fprintf(1,'a=%.6f m=%.6f s=%.6f\n',a,m,s);

Nt = randn(nexp,n);
S = zeros(nexp,n+1);

for i=1:nexp,
  S(i,:) = filter(1,[1, -a], [S0, m+s*Nt(i,:)]);
end


