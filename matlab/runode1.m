function runode1
% function runode1

%
% $Id: runode1.m,v 1.5 2011/03/26 09:20:41 patrick Exp $
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

global nui k sign2 sigv2

nui=1;
sign2=4*pi*nui/10;
sigv2=4*pi*nui/10;

ks=linspace(0,5*nui,10)';
y0 = linspace(0,1,10);
dy0 = linspace(0,1,10);

for i=1:10,
	k = ks(i);
	yi = [y0(i)*conj(y0(i));dy0(i)*conj(dy0(i)); ...
		y0(i)*conj(dy0(i))+conj(y0(i))*dy0(i)];
	[t,y]=ode45('ode1',[0; 25],yi);

	ya = ode1a(t,yi);
	if k~=0
		yas = 1/(4*nui*k^2)*(sigv2+k^2*sign2)+nui/k^2*sign2;
	else
		yas = Inf;
	end

	subplot(10,2,2*i-1),
	plot(t,y(:,2),t,ya(:,2))
	set(gca,'xlim',[min(t) max(t)]);
	ylabel(sprintf('|<d\\phi(k=%.2g,t)>|^2',k));
	if i==10,
		xlabel('t');
	end

	subplot(10,2,2*i),
	plot(t,y(:,1),t,ya(:,1),[min(t) max(t)],[yas yas]);
	set(gca,'xlim',[min(t) max(t)]);
	ylabel(sprintf('|<\\phi(k=%.2g,t)>|^2',k));
	if i==10,
		xlabel('t');
	end
end
axes('Position',[0 0 1 1],'Visible','off')
text(0.5,1, 'ode1', ...
	'HorizontalAlignment','center', ...
	'VerticalAlignment','top','Units','normalized');


orient tall
print -dpsc ode1
