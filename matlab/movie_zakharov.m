function M=movie_zakharov(solver,step)
% function M=movie_zakharov(solver,step)
% To save as mpeg file run mpgwrite(M,M(1).colormap,filename);

%
% $Id: movie_zakharov.m,v 1.3 2011/03/26 09:20:41 patrick Exp $
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

close all

mnx=min(solver.xj(:));
mxx=max(solver.xj(:));

mnn=min(solver.njs(:));
mxn=max(solver.njs(:));

mnE=min(abs(solver.Ejs(:)));
mxE=max(abs(solver.Ejs(:)));

for j=1:step:length(solver.times)

	subplot(211),
	plot(solver.xj, solver.njs(:,j))
	set(gca,'xlim',[mnx mxx]);
	set(gca,'ylim',[mnn mxn]);

	subplot(212),
	plot(solver.xj, abs(solver.Ejs(:,j)))
	set(gca,'xlim',[mnx mxx]);
	set(gca,'ylim',[mnE mxE]);

	M(j)=getframe(gcf);
	fprintf(1,'Image %d (%dx%dx%d)\n', j, size(M(j).cdata));

end
movie(M,1,2);

% mpgwrite(M,colormap(jet),mpgfilename)
