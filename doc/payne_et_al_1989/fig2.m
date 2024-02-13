function fig2
% function fig2

%
% $Id: fig2.m,v 1.2 2007/05/31 17:08:33 patrick Exp $
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

t=[];
Ek2=abs(sd('../../hdf/payne_et_al_1989_fig2_600','re-Ek',{t,[]})+i* ...
       sd('../../hdf/payne_et_al_1989_fig2_600','im-Ek',{t,[]})).^2;
nk2=abs(sd('../../hdf/payne_et_al_1989_fig2_600','re-nk',{t,[]})+i* ...
       sd('../../hdf/payne_et_al_1989_fig2_600','im-nk',{t,[]})).^2;

if 1
subplot(321), imagesc(Ek2'), set(gca,'xlim',[-20 20],'ylim',[0 150]);
subplot(322), imagesc(nk2'), set(gca,'xlim',[-20 20],'ylim',[0 150]);
else
plot(Ek2',[3 2 1]);
%set(gca,'xlim',[-20 20],'ylim',[0 5.5e-3]);
set(gca,'xlim',[-20 20]);
set(gca,'yticklabel',get(gca,'ytick'));
plot(nk2',[3 2 2]);
%set(gca,'xlim',[-20 20],'ylim',[0 5.5e-3]);
set(gca,'xlim',[-20 20]);
set(gca,'yticklabel',get(gca,'ytick'));
end


Ek2=abs(sd('../../hdf/payne_et_al_1989_fig2_1200','re-Ek',{t,[]})+i* ...
       sd('../../hdf/payne_et_al_1989_fig2_1200','im-Ek',{t,[]})).^2;
nk2=abs(sd('../../hdf/payne_et_al_1989_fig2_1200','re-nk',{t,[]})+i* ...
       sd('../../hdf/payne_et_al_1989_fig2_1200','im-nk',{t,[]})).^2;

if 1
subplot(323), imagesc(Ek2'), set(gca,'xlim',[-20 20],'ylim',[0 150]);
subplot(324), imagesc(nk2'), set(gca,'xlim',[-20 20],'ylim',[0 150]);
else
plot(Ek2',[3 2 3]);
set(gca,'xlim',[-20 20],'ylim',[0 5.5e-3]);
set(gca,'yticklabel',get(gca,'ytick'));
plot(nk2',[3 2 4]);
set(gca,'xlim',[-20 20],'ylim',[0 5.5e-2]);
set(gca,'yticklabel',get(gca,'ytick'));
end

Ek2=abs(sd('../../hdf/payne_et_al_1989_fig2_2400','re-Ek',{t,[]})+i* ...
       sd('../../hdf/payne_et_al_1989_fig2_2400','im-Ek',{t,[]})).^2;
nk2=abs(sd('../../hdf/payne_et_al_1989_fig2_2400','re-nk',{t,[]})+i* ...
       sd('../../hdf/payne_et_al_1989_fig2_2400','im-nk',{t,[]})).^2;

if 1
subplot(325), imagesc(Ek2'), set(gca,'xlim',[-20 20],'ylim',[0 150]);
subplot(326), imagesc(nk2'), set(gca,'xlim',[-20 20],'ylim',[0 150]);
else
plot(Ek2',[3 2 5]);
set(gca,'xlim',[-20 20],'ylim',[0 5.5e-3]);
set(gca,'yticklabel',get(gca,'ytick'));
plot(nk2',[3 2 6]);
set(gca,'xlim',[-20 20],'ylim',[0 5.5e-2]);
set(gca,'yticklabel',get(gca,'ytick'));
end


