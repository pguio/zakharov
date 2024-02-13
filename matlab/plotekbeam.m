function plotekbeam(s)
% function plotekbeam(s)

%
% $Id: plotekbeam.m,v 1.6 2011/03/26 09:20:41 patrick Exp $
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


t=[0, 150, 300, 450, 600, 750, 900, 1050, 1350, 1800, 3000, 3750];
ti=fix(t/s.msave)+1;

xlim=[-1 1];
coef=0.4;
index=1;

if (isfield(s,'Fejs'))
	mnp={311,312,313};
else
	mnp={211,212};
end

k=s.k*s.lambdae/s.chi;
ii=find(abs(k)<coef);
Ek=[];
nk=[];
for i=ti,
	Ek=[Ek, s.Eks(:,i).^2];
	nk=[nk, s.nks(:,i).^2];
end

ylim1=[min(min(Ek(ii,:))), max(max(Ek(ii,:)))];
ylim2=[min(min(nk(ii,:))), max(max(nk(ii,:)))];

for i=ti,

	if length(mnp)==3
		subplot(mnp{end}), 
		plot(s.vj/abs(s.vb), s.Fejs(:,i))
		set(gca,'xlim',[sign(s.vb)-0.5 sign(s.vb)+0.5]);
		xlabel('v/v_b','verticalalignment','top');
		ylabel('f_e(v/v_b)','verticalalignment','bottom');
	end

	subplot(mnp{1})
	plot(k, Ek(:,index))
	set(gca,'xlim',0.5*coef*xlim);
	set(gca,'ylim',ylim1);
	xlabel('k\lambda_e','verticalalignment','top');
	ylabel('|E(k)|^2','verticalalignment','bottom');

	title(sprintf('time = %.2e s', s.times(i)*s.tau);

	subplot(mnp{2})
	plot(k, nk(:,index))
	set(gca,'xlim',coef*xlim);
	set(gca,'ylim',ylim2);
	xlabel('k\lambda_e','verticalalignment','top');
	ylabel('|n(k)|^2','verticalalignment','bottom');

	orient tall
	if index==1,
		print ekbeam
	else
		print -append ekbeam
	end
	index=index+1;

	pause

end

