function plot_zakharov(solver)
% function plot_zakharov(solver)

%
% $Id: plot_zakharov.m,v 1.13 2011/03/26 09:20:41 patrick Exp $
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

if ~isfield(solver,'Fej')
	mnp={211,212};
else
	mnp={311,312,313};
end
	
switch solver.displaytype,
	case 'x',
		x=solver.xj*solver.chi;
		y1=real(solver.nj)/solver.nu;
		y2=[real(solver.Ej), imag(solver.Ej), abs(solver.Ej)]*solver.epsilon;
		f1='n(x)';
		f2='E(x)';
	case 'k',
		x=solver.k*solver.lambdae/solver.chi;
		y1=abs(solver.nk/solver.nu).^2;%/solver.N^2;
		y2=abs(solver.Ek*solver.epsilon).^2;%/solver.N^2;
		f1='|n(k)|^2';
		f2='|E(k)|^2';
end

subplot(mnp{1})
plot(x,y1)
set(gca,'xlim', [min(x) max(x)]);
ylabel(sprintf('%s', f1));
title(sprintf('time=%.2e/%.2e',  ...
	solver.ctime*solver.tau, solver.maxiter*solver.h*solver.tau))

subplot(mnp{2}) 
plot(x,y2)
set(gca,'xlim', [min(x) max(x)]);
ylabel(sprintf('%s', f2));

if length(mnp)==2,
	drawnow;
	return
end

subplot(mnp{3}) 


if isfield(solver,'dvb_vb') & isfield(solver,'vb')
if 0
	plot(solver.k, solver.Fej)
	set(gca,'xlim',solver.wpe./solver.vb*solver.chi.* ...
		[1-sign(solver.vb)*solver.dvb_vb 1+sign(solver.vb)*solver.dvb_vb]);
	ylabel('F_e(k)')
else
	x=solver.vj/abs(solver.vb);
	ii=find(x>sign(solver.vb)-0.5 & x<sign(solver.vb)+0.5);
	x=x(ii);
	y=solver.Fej(ii);
	plot(x, y)
  set(gca,'xlim',[sign(solver.vb)-0.5 sign(solver.vb)+0.5]);
	ylabel('F_e(v/vb)')
end
else
	plot(solver.k, solver.Fej)
	set(gca, 'xlim',[min(solver.k) max(solver.k)]);
	ylabel('F_e(k)')
end

drawnow;

