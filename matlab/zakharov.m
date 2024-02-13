function solver=zakharov(f,varargin)
% function solver=zakharov(f,varargin)

%
% $Id: zakharov.m,v 1.24 2011/03/26 09:20:42 patrick Exp $
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

solver=feval(f);
solver=zakharov_par_parsing(solver,varargin{:});
solver=zakharov_init(solver);

if solver.display
	plot_zakharov(solver);
end
view_zakharov(solver);

if solver.save
	solver.times=0;
	switch solver.displaytype,
		case 'x',
			solver.njs=solver.nj/solver.nu;
			solver.dnjs=solver.dnj/solver.nu/solver.tau;
			solver.Ejs=solver.Ej*solver.epsilon;
		case 'k',
			solver.nks=solver.nk/solver.nu;%/solver.N;
			solver.dnks=solver.dnk/solver.nu/solver.tau;%/solver.N;
			solver.Eks=solver.Ek*solver.epsilon;%/solver.N;
	end
	if isfield(solver,'Fej'),
		solver.Fejs=solver.Fej;
	end;
end
for i=1:solver.maxiter,

	solver.citer=i;
	solver=zakharov_time_integrate(solver);

	if mod(i,100)==0,
		s=sprintf(' Iter %6.d/%6.d: integrated to %.2f/%.2f (%.2f%%)',i, ...
			solver.maxiter,solver.ctime,solver.maxtime, ...
			100*solver.ctime/solver.maxtime);
		fprintf(1,'%s%s', s, char(repmat(8,1,length(s))));
	end
		
	if solver.save & mod(i,solver.msave)==0,
		solver.times=[solver.times solver.ctime];
		switch solver.displaytype
			case 'x',
				solver.njs=[solver.njs solver.nj/solver.nu];
				solver.dnjs=[solver.dnjs solver.dnj/solver.nu/solver.tau];
				solver.Ejs=[solver.Ejs solver.Ej*solver.epsilon];
			case 'k',
				solver.nks=[solver.nks solver.nk/solver.nu];%/solver.N];
				solver.dnks=[solver.dnks solver.dnk/solver.nu/solver.tau];%/solver.N];
				solver.Eks=[solver.Eks solver.Ek*solver.epsilon];%/solver.N];
		end
		if isfield(solver,'Fej'),
			solver.Fejs=[solver.Fejs solver.Fej];
		end
	end

	if solver.display & mod(i,solver.mdisplay)==0,
		plot_zakharov(solver);
	end

	% averaged fields
	if solver.citer>=solver.av_start & solver.citer<=solver.av_end,
		solver.av_nk = solver.av_nk+abs(solver.nk/solver.nu).^2;%/solver.N^2;
		solver.av_Ek = solver.av_Ek+abs(solver.Ek*solver.epsilon).^2;%/solver.N^2;
		solver.n_av = solver.n_av+1;
	end

	% is spectra
	if solver.citer>=solver.is_start & solver.citer<=solver.is_end,
		if solver.is_nfft==-1 & rem(solver.citer-solver.is_start, ...
				solver.is_stride)==0,
			it = fix((solver.citer-solver.is_start)/solver.is_stride)+1;
			solver.snkt(:,it) = solver.nk(solver.is_ks)/solver.nu;%/solver.N;
			solver.sEkt(:,it) = solver.Ek(solver.is_ks)*solver.epsilon;%/solver.N;
		end
		if solver.is_nfft>0 & rem(solver.citer-solver.is_start, ...
				solver.is_stride)==0,
			it = rem(fix((solver.citer-solver.is_start)/solver.is_stride), ...
				solver.is_nfft)+1;
			solver.snkt(:,it) = solver.nk(solver.is_ks)/solver.nu;%/solver.N;
			solver.sEkt(:,it) = solver.Ek(solver.is_ks)*solver.epsilon ...
				.*solver.k(solver.is_ks)/solver.chi;%/solver.N;
			if it==solver.is_nfft,
				solver.snkw = solver.snkw+(abs(fft(conj(solver.snkt)')).^2)';
				solver.sEkw = solver.sEkw+(abs(fft(conj(solver.sEkt)')).^2)';
				solver.n_is = solver.n_is+1;
				if 0
					subplot(211), plot(solver.is_fs,solver.snkw)
					subplot(212), plot(solver.is_fs,solver.sEkw)
					keyboard
				end
			end
		end
	end
end

if solver.maxiter,
	fprintf(1,' Iter %6.d/%6.d: integrated to %.2f/%.2f (%.2f%%)\n',i, ...
		solver.maxiter,solver.ctime,solver.maxtime, ...
		100*solver.ctime/solver.maxtime);
end

if solver.n_av>0,
	solver.av_nk = solver.av_nk/solver.n_av;
	solver.av_Ek = solver.av_Ek/solver.n_av;
end

if solver.n_is>0,
	for ik=1:length(solver.is_ks);
		solver.snkw(ik,:) = fftshift(solver.snkw(ik,:));
		solver.sEkw(ik,:) = fftshift(solver.sEkw(ik,:));
	end
	solver.snkw = solver.snkw/solver.n_is;
	solver.sEkw = solver.sEkw/solver.n_is;
end

if solver.display,
	plot_all_zakharov(solver);
end


