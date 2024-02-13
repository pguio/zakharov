function tst(type)
% function tst(type)

%
% $Id: tst.m,v 1.3 2011/03/26 09:20:42 patrick Exp $
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

if ~exist('type','var'),
	type='imagesc';
end

k=linspace(-20,20,1024)';
L=2*pi/(k(2)-k(1));
x=linspace(0,L,1024)';

fs=linspace(-0.5,0.5,165);
dt=1/(fs(2)-fs(1));
t=linspace(0,dt*164,165);


[X,T]=meshgrid(x,t);
y=cos(0.5*X-2*pi*1e-4*T);
y=cos(0.5*X-2*pi*1e-4*T)+cos(-4*X-2*pi*8e-4*T);
tm=mean(t);
y=exp(-(T-1/2*tm).^2/tm^2).*cos(0.5*X-2*pi*1e-4*T)+...
	exp(-(T-3/2*tm).^2/tm^2).*cos(-4*X-2*pi*8e-4*T);


subplot(311), 
imagesc(x,t,y)
axis xy
colorbar

subplot(312), 
% direct space transform
Y=zeros(size(y));
for i=1:size(y,1),
	Y(i,:)=fftshift(fftw(y(i,:),-1));
end
switch type,
	case 'imagesc',
		imagesc(k,t,abs(Y).^2), axis xy, colorbar
	case 'plot',
		plot(k,[mean(abs(Y).^2)-std(abs(Y)).^2; mean(abs(Y).^2)+std(abs(Y)).^2])
end

subplot(313)
i0=find(fs<0);
y1=zeros(size(y));
for i=1:size(y,2),
	% direct time transform
	Y=fftshift(fftw(y(:,i),1));
	% remove the negative frequencies
	Y(i0)=0.0; 
	% inverse time transform
	y1(:,i)=fftw(fftshift(Y),-1)/prod(size(Y));
end
% direct space transform
Y=zeros(size(y));
for i=1:size(y,1),
	Y(i,:)=fftshift(fftw(y1(i,:),-1));
end
switch type,
	case 'imagesc',
		imagesc(k,t,abs(Y).^2), axis xy, colorbar
	case 'plot',
		plot(k,[mean(abs(Y).^2)-std(abs(Y)).^2; mean(abs(Y).^2)+std(abs(Y)).^2])
end

