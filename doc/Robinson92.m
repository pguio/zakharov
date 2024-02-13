function Robinson92(figs)

opts = odeset('AbsTol',1e-8,'RelTol',1e-8);
seed = 11;
RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed));

if ~exist('figs') | isempty(figs)
  figs = 1:9;
end

for f=figs,

close all

switch f

case 1
Ga = 1; gb = 2; alpha = 1;
[t,y]=ode45(@twowaves,[0 40],[1 10],opts,Ga,gb,alpha);
plot(t,log10(y),t,log10(sum(y')));
title('Figure 1')

case 2
Ga = 1; gb = 2; gc = 10; alpha = 1;
[t,y]=ode45(@threewaves,[0 40],[1 10 1],opts,Ga,gb,gc,alpha);
plot(t,log10(y));
title('Figure 2a')

figure
Ga = 1; gb = 0.3; gc = 0.6; alpha = 1;
[t,y]=ode45(@threewaves,[0 20],[1 10^.3 10^.6],opts,Ga,gb,gc,alpha);
plot(t,log10(y));
title('Figure 2b')

case 3
Ga = 0; gb = 1; alpha = 1;
[t,y]=ode15s(@twowaves,[0 20],[10 10^-5],opts,Ga,gb,alpha);
plot(t,log10(y));
title('Figure 3')

case 4
Ga = 0; gb = 1; gc = 20; alpha = 1;
[t,y]=ode15s(@threewaves,[0 20],[10 10^-5 10^-5],opts,Ga,gb,gc,alpha);
plot(t,log10(y));
title('Figure 4')

case 5
Ga = 0; gb = 1; gc = 1; alpha = 1;
[t,y]=ode15s(@threewaves,[0 20],[10^0.6 10^-5 10^-5],opts,Ga,gb,gc,alpha);
plot(t,log10(y));
title('Figure 5a')

figure
Ga = 0; gb = 1; gc = 1; alpha = 1;
[t,y]=ode15s(@threewaves,[0 20],[10^1.25 10^-5 10^-5],opts,Ga,gb,gc,alpha);
plot(t,log10(y));
title('Figure 5b')

case 6
dt = 0.1; T = 1000;
Ga = 1; dGa = 1;
Gas = Ga+dGa*(2*rand(fix(T/dt)+1,1)-1); 
gb = 1;
alpha = 1;
S = 0;
[t,y]=ode45(@twowavesstoch,[0 T],[1 1],opts,Gas,gb,alpha,dt,S);
plot(1e-2*t,log10(y(:,1)));
title('Figure 6a')

figure
plot(1e-2*t,Vfun(t,y,Ga,gb,alpha,dt));
title('Figure 6b')

case 7
dt = 0.1; T = 5000;
Ga = 1; dGa = 0.02;
Gas = Ga+dGa*(2*rand(fix(T/dt)+1,1)-1); 
gb = 1;
alpha = 1;
S = 0;
[t,y]=ode45(@twowavesstoch,[0 T],[1 1],opts,Gas,gb,alpha,dt,S);
plot(1e-2*t,log10(y(:,1:2)));
title('Figure 7a')

figure
Ga = 1; dGa = 0.02;
Gas = Ga+dGa*(2*rand(fix(T/dt)+1,1)-1);
gb = 1;
alpha = 1;
S = 0.02;
[t,y]=ode45(@twowavesstoch,[0 T],[1 1],opts,Gas,gb,alpha,dt,S);
plot(1e-2*t,log10(y(:,1:2)));
title('Figure 7b')

case 8
dt = 0.1; T = 100;
Ga = 0.1; dGa = 2;
Gas = Ga+dGa*(2*rand(fix(T/dt)+1,1)-1);
gb = 1;
gc = 10;
alpha = 1;
S = 0;
[t,y]=ode45(@threewavesstoch,[0 T],[1 1e-2 1e-3],opts,Gas,gb,gc,alpha,dt,S);
plot(t,log10(y(:,1:3)));
title('Figure 8');

case 9
dt = 0.1; T = 100;
Ga = 0.1; dGa = 2;
Gas = Ga+dGa*(2*rand(fix(T/dt)+1,1)-1);
gb = 1;
gc = 1;
alpha = 1;
S = 0;
[t,y]=ode45(@threewavesstoch,[0 T],[1 1e-2 1e-2],opts,Gas,gb,gc,alpha,dt,S);
plot(t,log10(y(:,1:3)));
title('Figure 9');

end
if length(figs)>1,
pause
end
end

function yp = twowaves(t,y,Ga,gb,alpha)

yp = zeros(size(y));
yp(1) = Ga*y(1)-alpha*y(1)*y(2);
yp(2) = -gb*y(2)+alpha*y(1)*y(2);


function yp = threewaves(t,y,Ga,gb,gc,alpha)

yp = zeros(size(y));
yp(1) = Ga*y(1)-alpha*(y(1)*y(2)-y(2)*y(3)+y(1)*y(3));
yp(2) = -gb*y(2)+alpha*(y(1)*y(2)-y(2)*y(3)+y(1)*y(3));
yp(3) = -gc*y(3)+alpha*(y(1)*y(2)-y(2)*y(3)+y(1)*y(3));


function yp = twowavesstoch(t,y,Gas,gb,alpha,dt,S)

Ga = Gas(fix(t/dt)+1);

yp = zeros(size(y));
yp(1) = Ga*y(1)-alpha*y(1)*y(2)+S;
yp(2) = -gb*y(2)+alpha*y(1)*y(2)+S;

function yp = threewavesstoch(t,y,Gas,gb,gc,alpha,dt,S)

Ga = Gas(fix(t/dt)+1);

yp = zeros(size(y));
yp(1) = Ga*y(1)-alpha*(y(1)*y(2)-y(2)*y(3)+y(1)*y(3))+S;
yp(2) = -gb*y(2)+alpha*(y(1)*y(2)-y(2)*y(3)+y(1)*y(3))+S;
yp(3) = -gc*y(3)+alpha*(y(1)*y(2)-y(2)*y(3)+y(1)*y(3))+S;

function V = Vfun(t,y,Ga,gb,alpha,dt)

V = alpha*(y(:,1)+y(:,2))-gb*log(y(:,1))-Ga.*log(y(:,2));
