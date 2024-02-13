function [kn, win]=damping

L = 100;
N = 4096;

ne = 5e11;
Te = 3000;
nuec = 100;

mi = 16;
Ti = 1000;

nb_ne = 3e-5;
ub = 4.5e6;
dub_ub = 0.5;

%%%%%%%

k = 2*pi/L*[-N/2:-1, 1:N/2-1];
ip = find(k > 0);
in = find(k < 0);

si=SIconstants;

me = si.me;
eta = (Te+3*Ti)/Te;
we = sqrt(ne*si.e^2/(si.eps0*me));
ve = sqrt(si.kb*Te/me);
lambdae = sqrt(si.eps0*si.kb*Te/(ne*si.e^2));

mi = mi*si.amu;

dub = ub*dub_ub;

%%%%%%%%

alpha = k*lambdae;
beta = sqrt(1+3*alpha.^2);

wr = we*beta;

plot((wr-we)/(2*pi), k)
set(gca,'ylim',[10 40]);
pause

wim = sqrt(pi/8)*we./abs(alpha).^3.*(1+3*alpha.^2).*exp(-0.5./alpha.^2-1.5);

wib = sqrt(pi/8)*we./alpha.^2*nb_ne*(ve/dub)^2.*beta.*abs(k)./k.* ...
	(beta./alpha*(ve/dub)-ub/dub).*exp(-0.5*(beta./alpha*(ve/dub)-ub/dub).^2);

wi = 0.5*nuec+wim+wib;

fm = 1/sqrt(2*pi)/ve*exp(-0.5./alpha.^2-1.5);
fb = nb_ne/sqrt(2*pi)/dub*exp(-0.5*((wr./k-ub)/dub).^2);
f = fm+fb;

if 1
	dkdv = -k.^3/we^2.*wr./abs(k);
	df=zeros(size(f));
	df(ip) = derive(k(ip), f(ip)).*dkdv(ip);
	df(in) = derive(k(in), f(in)).*dkdv(in);
	wi2 = -pi/2*(we./k).^2.*wr.*df;
else
	df(ip) = derive(k(ip), f(ip));
	df(in) = derive(k(in), f(in));
	wi2 = pi/2*sign(k)*we^2.*(1+3*alpha.^2).*df;
end
wi2 = 0.5*nuec+wi2;

%subplot(111), plot(alpha, wr), pause
subplot(211), plot(alpha, wi, alpha, wi2)
subplot(212), plot(alpha, wi, alpha, wi2)
set(gca,'xlim',[-0.2 0.2]);
set(gca,'ylim',[0 500]);

pause

%%%%%%%%

tau = 3/2*(mi/(me*eta))*(1/we);
chi = 3/2*sqrt(mi/(me*eta))*lambdae;

%%%%%%

kn = 3/2*sqrt(mi/(me*eta))*lambdae*k;

wimn = sqrt(pi/8)*(3/2)^4*(mi/(me*eta))^2.5./abs(kn).^3.* ...
	(1+4/3*kn.^2*me*eta/mi).*exp(-9/8*mi/(me*eta)./kn.^2-3/2);

wibn = sqrt(pi/8)*(3/2)^3*(mi/(me*eta))^2./kn.^2.*sign(kn)* ...
	nb_ne*(ve/dub)^2.*beta.*(beta./alpha*(ve/dub)-ub/dub).* ...
	exp(-0.5*(beta./alpha*(ve/dub)-ub/dub).^2);

win = 0.5*nuec*tau+wimn+wibn;

wimn1 = sqrt(pi/8)*(3/2)^4*(mi/(me*eta))^2.5./abs(kn).^3.* ...
	exp(-(9/8)*mi/(me*eta)./kn.^2-3/2);

x = we/ub./k;
wibn1 = -sqrt(pi/8)*(3/2)*mi/(me*eta)*(x/dub_ub).^2.*sign(x)* ...
	nb_ne/dub_ub.*(1-x).*exp(-0.5*((1.0-x)/dub_ub).^2);

win1 = 0.5*nuec*tau+wimn1+wibn1;

subplot(211), plot(alpha, wi*tau, alpha, win, alpha, win1)
subplot(212), plot(alpha, wi*tau, alpha, win, alpha, win1)
set(gca,'xlim',[-0.2 0.2]);
set(gca,'ylim',[-0.5 1]);
pause

subplot(111), plot(k, wi*tau, k, win)
%subplot(111), plot(k, wi*tau-win)


function dy=derive(x,y)
dy = zeros(size(y));
dy(1) = (y(2)-y(1))/(x(2)-x(1));
dy(2:end-1) = (y(3:end)-y(1:end-2))./(x(3:end)-x(1:end-2));
dy(end) = (y(end-1)-y(end))/(x(end-1)-x(end));

