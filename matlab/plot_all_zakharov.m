function plot_all_zakharov(solver,wtk)
% function plot_all_zakharov(solver,wtk)

%
% $Id: plot_all_zakharov.m,v 1.7 2011/03/26 09:20:41 patrick Exp $
%
% Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
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

if ~exist('wtk','var') | isempty('wtk')
	wtk=0;
end

switch (solver.displaytype),
	case 'x',
		x=solver.xj;         tx='x';
		n=real(solver.njs)'; tn='n(x,\tau)';
		E=abs(solver.Ejs)';  tE='|E(x,\tau)|';
	case 'k';
		x=solver.k;            tx='k';
		switch wtk,
			case 0,
				n=abs(solver.nks').^2; tn='|n(k,\tau)|^2';
				E=abs(solver.Eks').^2; tE='|E(k,\tau)|^2';
			otherwise,
				nks=wkt(solver.nks,solver);
				n=real(nks').^2; tn='|n(k,\tau)|^2';
				E=abs(solver.Eks').^2; tE='|E(k,\tau)|^2';
		end
end
y=solver.times; ty='\tau';

subplot(211), 
if 1
imagesc(x, y, n);
axis xy
colorbar
ylabel(ty)
else
plot(x,n(end,:))
end
title(tn)

subplot(212), 
if 1
imagesc(x, y, E);
axis xy
colorbar
ylabel(ty)
else
plot(x,E(end,:))
end
xlabel(tx)
title(tE)

function [nks]=wkt(nks, solver)

njs=zeros(size(nks));

% inverse space transform
for it=1:size(nks,2),
	njs(:,it)=IFFT(nks(:,it),solver.M,solver.N);
	%subplot(111), plot([real(njs(:,it)), imag(njs(:,it))]); drawnow
end

nks=zeros(size(njs));

dt=solver.times(2)-solver.times(1);
fs=1/dt*[-length(solver.times)/2:length(solver.times)/2-1];

i0=find(fs<0);
njs1=zeros(size(njs));
for ik=1:size(njs,1),
	% direct time transform
	nk=fftshift(fftw(njs(ik,:),1));
	% remove the negative frequencies
	nk(i0)=0.0;
	%subplot(111), plot(fs,abs(nk)); drawnow
	% inverse time transform
	njs1(ik,:)=fftw(fftshift(nk),-1)/prod(size(nk));
end
% direct space transform
for it=1:size(njs,2),
	nks(:,it)=FFT(njs1(:,it),solver.M,solver.N);
end


