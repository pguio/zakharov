function cmplw
% function cmplw

%
% $Id: cmplw.m,v 1.5 2011/03/26 09:20:41 patrick Exp $
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

if 0

	dt = [];
	ik2 = [];
	uk2 = [];
	dk2 = [];
	%for i=[1:10 20 30 40]
	for i=[1:10]
		eval(['load s' num2str(i)]);
		dt = [dt s.H];
		ik2 = [ik2; s.ispec.sk2(1:20:end)];
		uk2 = [uk2; s.uspec.sk2(1:20:end)];
		dk2 = [dk2; s.dspec.sk2(1:20:end)];
	end
	subplot(311), plot(dt, ik2)
	subplot(312), plot(dt, uk2)
	subplot(313), plot(dt, dk2)

else

	dt = [];
	ik2 = [];
	uk2 = [];
	dk2 = [];
	j = 1;
	%for i=[1:10 20 30 40]
	for i=[1:10]
		eval(['load s' num2str(i)]);
		dt = [dt s.H];
		leg{j} = sprintf('h=%.2e', s.H); %num2str(s.H,3);
		ik2 = [ik2; s.ispec.sk2];
		uk2 = [uk2; s.uspec.sk2];
		dk2 = [dk2; s.dspec.sk2];
		subplot(13,2,2*j-1), 
		plot(s.ispec.k, s.ispec.sk2)
		set(gca,'ylim',[0 6e10]);
		set(gca,'xlim',[10 40]);
		subplot(13,2,2*j), 
		plot(s.uspec.k, s.uspec.sk2, abs(s.dspec.k), s.dspec.sk2)
		set(gca,'ylim',[0 6e7]);
		set(gca,'xlim',[10 40]);
		j = j+1;
	end
	subplot(111), plot(s.ispec.k, ik2, 'linewidth', 1.5)
	xlabel('k [m^{-1}]')
	ylabel('<|\delta n(k)|^2>')
	set(gca,'xlim',[10 40]);
	legend(leg,1)
	%subplot(312), plot(s.uspec.k, uk2)
	%subplot(313), plot(abs(s.dspec.k), dk2)

end
