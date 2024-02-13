function tfanalysis
% function tfanalysis

%
% $Id: tfanalysis.m,v 1.4 2011/03/26 09:20:42 patrick Exp $
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

close all

if ~exist('nkt.mat','file')
	load nkt.out
	t=nkt(:,1);
	nkt=nkt(:,2:2:end)+i*nkt(:,3:2:end);
	save nkt nkt t
else
	load nkt
end

sig=nkt(:,1);

for i=1:20,
	[tfr,t,f]=tfrsp(sig([(i-1)*1024+1:i*1024]));
	subplot(5,4,i), imagesc(t,f,tfr);
end

if ~exist('Ekt.mat','file')
	load Ekt.out
	t=Ekt(:,1);
	Ekt=Ekt(:,2:2:end)+i*Ekt(:,3:2:end);
	save Ekt Ekt t
else
	load Ekt
end


sig=Ekt(:,2);

[tfr,t,f]= tfrstft(Ekt(1:1024,7));  
if 0
f=fftshift(f);
for i=1:size(tfr,2),
	tfr(:,i)=fftshift(tfr(:,i));
end
imagesc(t, f, abs(tfr).^2)
else
imagesc(t, f(1:end/2), abs(tfr(1:end/2,:)).^2)
end
axis xy
xlabel('Time [s]');
ylabel('Frequency [Hz]');

