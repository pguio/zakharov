function test_epsilons
% function test_epsilons

%
% $Id: test_epsilons.m,v 1.1 2001/02/07 13:05:38 patricg Exp $
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

T = 20/(2*pi/10);
h=[T/400 T/375 T/350 T/325 T/300 T/275 T/250 T/225 T/200];
it=[300:-10:200];
h=T./it;

epsilon = zeros(size(h));

for i=1:length(h),
	solver{i}=zakharov('A','N',64, ...
		'display',0,'save',1,'msave',10,'h',h(i));
	[h(i),epsilon(i)] = zakharov_epsilon(solver{i});
	[t, epsilons] = zakharov_epsilons(solver{i});
	ts{i}=t;
	es{i}=epsilons;
end

subplot(211)
plot(h,epsilon,'-o')
set(gca,'xlim',[min(h) max(h)],'ylim',[min(epsilon) max(epsilon)])
xlabel('h')
title('\epsilon_h(D)')

subplot(212)
for i=1:length(h)
	plot(ts{i},es{i},'-o')
	hold on
end
set(gca,'xlim',[0 T])
xlabel('t');
title('\epsilon_h(t)')

