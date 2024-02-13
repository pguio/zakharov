function fig2
% function fig2

%
% $Id: fig2.m,v 1.3 2007/05/27 14:40:33 patrick Exp $
%
% Copyright (c) 2000 Patrick Guio <patrick@phys.uit.no>
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

close all

if exist('fig2.mat')
  load fig2
else
	MA = [5 20 40];
	epsilonA = zeros(size(MA));
	i = 1;
	for M=MA,
		solverA=zakharov('A','N',128,'display',0,'maxiter',0,'h',0);
		solverA.Ej = IFFT(FFT(solverA.Ej,solverA.N,M),solverA.N,M);
		[h, epsilonA(i)] = zakharov_epsilon(solverA);
		i=i+1;
	end

	MB = [10 20 40 80 170 340];
	epsilonB = zeros(size(MB));
	i = 1;
	for M=MB,
		solverB=zakharov('B','N',1024,'display',0,'maxiter',0,'h',0);
		solverB.Ej = IFFT(FFT(solverB.Ej,solverB.N,M),solverB.N,M);
		[h, epsilonB(i)] = zakharov_epsilon(solverB);
		i=i+1;
	end

	save fig2 epsilonA epsilonB MA MB
end

epsilonA = log10(epsilonA);
epsilonB = log10(epsilonB);
h=plot(MA,epsilonA,'-o',MB,epsilonB,'-d');
%h=semilogy(MA,epsilonA,'-o',MB,epsilonB,'-d');
set(h(:),'MarkerSize',5)
set(gca,'xlim',[0 360],'ylim',log10([1e-10 1e0]));
xlabel('NUMBER OF MODES M');
ylabel('LOG_{10} \epsilon');
legend('E_{max}=1.0','E_{max}=10.0')

orient landscape; set(gcf,'PaperOrientation','portrait');
if exist('exportfig'),
	exportfig(gcf,'fig2');
else
	print -deps fig2
end

