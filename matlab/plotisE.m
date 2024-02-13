function [f,is,k,int_is]=plotisE(s,n_av,plotcmd,mode)
% function [f,is,k,int_is]=plotisE(s,n_av,plotcmd,mode)

%
% $Id: plotisE.m,v 1.18 2011/03/26 09:20:41 patrick Exp $
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

if ~isfield(s,'sEkw') 

  if exist('n_av')~=1,
		n_av = 1;
	end

	if n_av~=1 & rem(size(s.sEkt,2), n_av)~=0
		error(sprintf('n_av must be a multiple of %d',size(s.sEkt,2)))
	end

	if ~exist('mode'),
		mode=0;
	end

	if mode==0,
		is=zeros(size(s.sEkt,1), size(s.sEkt,2)/n_av);
	else
		is=zeros(size(s.sEkt,1), size(s.sEkt,2)/n_av, n_av);
	end	
	for ik=1:length(s.is_ks),
		it=[1:size(s.sEkt,2)/n_av];
		specs=zeros(size(s.sEkt,2)/n_av,n_av);
		for k=1:n_av,
			% fprintf(1,'min(it)=%d max(it)=%d\n',min(it),max(it));
			specs(:,k) = fftshift(abs(fft(conj(s.sEkt(ik,it))')).^2);
			if mode==1,
				is(ik,:,k) = specs(:,k);
			end
			it = it + size(s.sEkt,2)/n_av;
		end
		if n_av==1,
			is(ik,:) = specs';
		else
			if mode==0,
				is(ik,:) = mean(specs');
			end
		end
		is(ik,:) = is(ik,:)*s.k(s.is_ks(ik))^2;
	end

	f=1/s.h*[-size(is,2)/2:size(is,2)/2-1]/size(is,2);
	is = is*(s.epsilon*s.tau)^2; % [(v m^-1 s)^2]

else

	is = s.sEkw;
	f = s.is_fs;

end


n=s.n_ks;
n2=fix(s.n_ks/2);

if ~exist('plotcmd')
	plotcmd='plot';
end

int_is=zeros(size(n));
k=s.k(s.is_ks)/s.chi;

for i=1:n2,
	subplot(n2,2,2*i-1)
	if s.k(s.is_ks(i))<0.0,
		if mode==0,
			feval(plotcmd,-f,is(i,:),'-')
		else
			imagesc(-f,1:n_av,squeeze(is(i,:,:))');
		end
	else
		if mode==0,
			feval(plotcmd,f,is(i,:),'-')
		else
			imagesc(f,1:n_av,squeeze(is(i,:,:))');
		end
	end

	x=1e3*f(:);
	if mode==0,
		y=is(i,:);
		int_is(i)=integrate(x(:), y(:));
	end

	if mode==0,
		set(gca,'ylim',[min(is(i,:)) max(is(i,:))]);
	else
		set(gca,'clim',[min(is(i,:)) max(is(i,:))]);
	end
	%set(gca,'xlim',0.5*[min(f) max(f)]);
	title(sprintf('k=%.1f m^{-1}',s.k(s.is_ks(i))/s.chi));
	if i~=n2,
		set(gca,'xticklabel',[])
	end
	xlabel('Frq [kHz]')

	subplot(n2,2,2*i)
	if s.k(s.is_ks(n2+i))<0.0,
		if mode==0,
			feval(plotcmd,-f,is(n2+i,:),'-')
		else
			imagesc(-f,1:n_av,squeeze(is(n2+i,:,:))');
		end
	else
		if mode==0,
			feval(plotcmd,f,is(n2+i,:),'-')
		else
			imagesc(f,1:n_av,squeeze(is(n2+i,:,:))');
			%plot(squeeze(is(n2+i,:,:)))
		end
	end

	x=1e3*f(:);
	if mode==0,
		y=is(n2+i,:);
		int_is(n2+i)=integrate(x(:), y(:));
	end

	if mode==0,
		set(gca,'ylim',[min(is(n2+i,:)) max(is(n2+i,:))]);
	else
		set(gca,'clim',[min(is(n2+i,:)) max(is(n2+i,:))]);
	end
	%set(gca,'xlim',0.5*[min(f) max(f)]);
	title(sprintf('k=%.1f m^{-1}',s.k(s.is_ks(n2+i))/s.chi));
	if i~=n2,
		set(gca,'xticklabel',[])
	end
	if i==n2,
		xlabel('Frq [kHz]')
	end
end
suptitle('<k^2|E(k,f)|^2>');
