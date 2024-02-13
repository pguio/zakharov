% function solver=zakharov_sources(solver)
function solver=zakharov_sources(solver)

%
% $Id: zakharov_sources.m,v 1.13 2011/03/26 09:20:42 patrick Exp $
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

if isfield(solver,'delta_n') | isfield(solver,'delta_E'),

	nui = solver.nui;
	if 1
		nue = solver.nueL;
	else
		nue = solver.nue;
	end
	h = solver.h;
	k = solver.k;

	ct = solver.ct;
	st = solver.st;

% Calculate new density fluctuation with random phase
	delta_nhk = solver.delta_n.*exp(i*2*pi*rand(size(k)));
% Make the density spectral fluctuation  such that
% delta_nhk(-k)=conj(delta_nhk(k))
	ip = find(solver.signed_index>0);
	im = find(solver.signed_index<0);
	delta_nhk(ip) = flipud(conj(delta_nhk(im(2:end))));
	ii = find(solver.signed_index==0);
	delta_nhk(ii) = 0.5*(delta_nhk(ii)+conj(delta_nhk(ii)));

	if 0
		subplot(211), plot(real(IFFT(delta_nhk))), title('delta_nhk'),
		subplot(212), plot(imag(IFFT(delta_nhk))), pause
	end

% Calculate new electric field fluctuation with random phase
	delta_Ehk=solver.delta_E.*exp(i*2*pi*rand(size(k)));

%	solver.snk = 2*nui.*solver.delta_nk.*exp(-nui*h).*st - ...
%		2*nui*(h/2).*(solver.delta_nk.*exp(-nui*h).*(nui.*st-ct)-delta_nhk);
	solver.snk = delta_nhk;

	if 0
		subplot(211), plot(real(IFFT(solver.snk))), title('solver.snk'),
		subplot(212), plot(imag(IFFT(solver.snk))), pause
	end

%	solver.sEk = k.*nue.*(solver.delta_Ek.*exp(-(i*k.^2+nue)*h)+delta_Ehk);
	solver.sEk = delta_Ehk;

%	solver.sdnk = 2*nui.*(delta_nhk-solver.delta_nk.*exp(-nui*h))- ...
%		4*nui.^2*(h/2).*(solver.delta_nk.*exp(-nui*h)+delta_nhk);
	solver.sdnk = zeros(size(k));

	if 0
		subplot(211), plot(real(IFFT(solver.sdnk))), title('solver.sdnk'),
		subplot(212), plot(imag(IFFT(solver.sdnk))), pause
	end

	solver.delta_nk = delta_nhk;
	solver.delta_Ek = delta_Ehk;

	solver.snk(solver.ikeq0) = 0;
	solver.sEk(solver.ikeq0) = 0;
	solver.sdnk(solver.ikeq0) = 0;

	if 0,
		subplot(311), plot(solver.k,real(solver.snk),solver.k,imag(solver.snk));
		subplot(312), plot(solver.k,real(solver.sEk),solver.k,imag(solver.sEk));
		subplot(313), plot(solver.k,real(solver.sdnk),solver.k,imag(solver.sdnk));
		keyboard
		pause
	end

else

	solver.snk = 0;
	solver.sEk = 0;
	solver.sdnk = 0;

end
