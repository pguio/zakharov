function teststint(N,N1,seeds,nu,sigma)
% function teststint(N,N1,seeds,nu,sigma)
% try to run
% teststint(500,5,7+[1:10],50,10) 

%
% $Id: teststint.m,v 1.3 2011/03/26 09:20:42 patrick Exp $
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

T=1;
dt = T/N;
t = [0,dt:dt:1];

dt1 = dt/N1;
t1=[0,dt1:dt1:dt];

NS = length(seeds);

toto = sigma*ones(NS,N+1);
ito = sigma*ones(NS,N+1);
strat = sigma*ones(NS,N+1);

ito1 = sigma*ones(NS,N+1);
strat1 = sigma*ones(NS,N+1);

for s = 1:NS,

	randn('state',seeds(s))

	dB = sqrt(dt)*randn(N1,N);

	for i=2:N
		ito(s,i+1) = ito(s,i)*exp(-nu*dt)+ ...
			exp(-nu*dt)*dB(1,i)*sigma*sqrt(2*nu);
		strat(s,i+1) = strat(s,i)*exp(-nu*dt)+ ...
			exp(-nu*dt/2)*dB(1,i)*sigma*sqrt(2*nu);
		toto(s,i+1) = toto(s,i)*exp(-nu*dt)+dB(1,i)*sigma*sqrt(2*nu);

		ito1(s,i+1) = ito1(s,i)*exp(-nu*dt)+ ...
			exp(-nu*dt)*sum(exp(nu*t1(1:end-1)).*dB(:,i)'/sqrt(N1))*sigma*sqrt(2*nu);
		strat1(s,i+1) = strat1(s,i)*exp(-nu*dt)+ ...
			exp(-nu*dt)*sum(exp(nu*(t1(1:end-1)+t1(2:end)/2)).*dB(:,i)'/sqrt(N1))* ...
			sigma*sqrt(2*nu);
	end
	
end


fprintf(1,'nu*dt=%3g\n', nu*dt);
fprintf(1,'mean(toto^2)=%3g\n', mean(mean(toto.^2)));
fprintf(1,'mean(ito^2)=%3g\n', mean(mean(ito.^2)));
fprintf(1,'mean(strat^2)=%3g\n', mean(mean(strat.^2)));
fprintf(1,'mean(ito1^2)=%3g\n', mean(mean(ito1.^2)));
fprintf(1,'mean(strat1^2)=%3g\n', mean(mean(strat1.^2)));

subplot(211), plot(t, mean(ito.^2), t, mean(strat.^2), t, mean(toto.^2));
subplot(212), plot(t, mean(ito1.^2), t, mean(strat1.^2));



