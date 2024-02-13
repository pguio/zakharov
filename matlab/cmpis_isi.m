function cmpis_isi
% function cmpis_isi

%
% $Id: cmpis_isi.m,v 1.3 2011/03/26 09:20:41 patrick Exp $
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

fi=sd('hdf/is4096i','dkf2');
imagesc(fi', [1 3 1])  
axis ij
fi=sd('hdf/is4096i','nkf2');
imagesc(fi', [1 3 2])  
fi=sd('hdf/is4096i','ukf2');
imagesc(fi', [1 3 3])  

return

ii=1:50:501;

figure
fi=sd('hdf/is4096i','nkf2',{[],[]});
plot(squeeze(subvolume(fi,{ii,[]}))',[2 1 1])  
f=sd('hdf/is4096','nkf2',{[],[],[]});
plot(squeeze(subvolume(mean(f),{ii,[]})'),[2 1 2])

figure
fi=sd('hdf/is4096i','ukf2',{[],[]});
plot(squeeze(subvolume(fi,{ii,[]}))',[2 1 1])  
f=sd('hdf/is4096','ukf2',{[],[],[]});
plot(squeeze(subvolume(mean(f),{ii,[]})'),[2 1 2])

figure
fi=sd('hdf/is4096i','dkf2',{[],[]});
plot(squeeze(subvolume(fi,{ii,[]}))',[2 1 1])  
f=sd('hdf/is4096','dkf2',{[],[],[]});
plot(squeeze(subvolume(mean(f),{ii,[]})'),[2 1 2])

