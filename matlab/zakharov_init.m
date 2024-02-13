function solver=zakharov_init(solver)
% function solver=zakharov_init(solver)

%
% $Id: zakharov_init.m,v 1.32 2011/03/26 09:20:42 patrick Exp $
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

% Initialise seed
rand('state',0);

% Initialise to 1 coefficient for normalised/physical variables
% if not defined
if ~isfield(solver,'tau'), solver.tau=1; end
if ~isfield(solver,'chi'), solver.chi=1; end
if ~isfield(solver,'epsilon'), solver.epsilon=1; end
if ~isfield(solver,'nu'), solver.nu=1; end
if ~isfield(solver,'lambdae'), solver.lambdae=1; end
if ~isfield(solver,'type'), solver.type='efluct'; end

% Grid points xj in configuration space
solver.xj = [0:solver.N-1]' * solver.L / solver.N;
solver.xj = solver.xj - solver.L / 2.0;

% M < N/3 
if ~isfield(solver,'M') solver.M = floor(solver.N/3); end

% Calculate signed_index of k's
solver.signed_index=[-solver.N/2:solver.N/2-1]';
% Index of k's to set to zero (large wave-numbers)
solver.ikeq0=find(solver.signed_index>solver.M|solver.signed_index<-solver.M);

% Grid points k in fourier space
solver.k = 2 * pi * [-solver.N/2:solver.N/2-1]' / solver.L;

% Time initialisation
solver.ctime=0;

% Init E(k), E(x), n(k), n(x), dn(k), dn(x)
switch lower(solver.type),
	case 'solitons'
		solver=soliton_init(solver);
	case 'efluct'
		solver=efluct_init(solver);
	case 'egauss'
		solver=egauss_init(solver);
	case 'tfluct'
		solver=thermal_fluctuations_init(solver);
		solver=tfluct_init(solver);
	case 'efluct_grade'
		solver=efluct_grade_init(solver);
end

if ~isfield(solver,'D'), solver.D = solver.L; end

% Calculation of the cos and sin terms in Eq. 10 
k = solver.k;
nui = solver.nui;
h = solver.h;
frq=k.^2-nui.^2;
solver.ct=ones(size(k));
solver.st=repmat(h,size(k));
if 1
ip=find(frq>0);
solver.ct(ip) = cos(sqrt(frq(ip))*h);
solver.st(ip) = sin(sqrt(frq(ip))*h)./sqrt(frq(ip));
im=find(frq<0);
solver.ct(im) = cosh(sqrt(-frq(im))*h);
solver.st(im) = sinh(sqrt(-frq(im))*h)./sqrt(-frq(im));
%solver.ct = cos(sqrt(frq)*h);
%ii=find(frq~=0);
%solver.st(ii) = sin(sqrt(frq(ii))*h)./sqrt(frq(ii));
else
ip=find(frq>0);
solver.ct(ip) = cos(sqrt(frq(ip))*h);
solver.st(ip) = sin(sqrt(frq(ip))*h)./sqrt(frq(ip));
end

% Initialise variables for averaged fields
if ~isfield(solver,'av_start') | ~isfield(solver,'av_end')
	solver.av_start = -1;
	solver.av_end = -1;
end
solver.av_Ek = zeros(size(solver.Ek));
solver.av_nk = zeros(size(solver.nk));
solver.n_av = 0;

% Test if kmin, kmax and nbk exists and create the array is_ks
% (Wave vectors to consider to calculate IS spectra)
if isfield(solver,'is_kmin') & isfield(solver,'is_kmax') & ...
	isfield(solver,'is_nbk'),
	[diff,kmn]=min(abs(solver.k-solver.is_kmin*solver.chi));
	[diff,kmx]=min(abs(solver.k-solver.is_kmax*solver.chi));
	solver.is_ks=fix(linspace(kmn,kmx,solver.is_nbk));
end

% Test if the selected k's are within [-M, M]
if isfield(solver,'is_ks'),
	ii=find(solver.is_ks-solver.N/2>=solver.M);
	if ~isempty(ii),
		error(sprintf('Some of the IS k''s are outside the range [-M, M]\n %d',ii));
	end
end

% Initialise variables for spectra,
if ~isfield(solver,'is_start') | ~isfield(solver,'is_end') | ...
	~isfield(solver,'is_stride') | ~isfield(solver,'is_ks'),
	solver.is_start = -1;
	solver.is_end = -1;
else	
	solver.n_ks = length(solver.is_ks);
	if ~isfield(solver,'is_nfft'),
		solver.is_nfft=-1;
		solver.is_nt = fix((solver.is_end-solver.is_start)/solver.is_stride);
		solver.snkt=zeros(solver.n_ks,solver.is_nt);
		solver.sEkt=zeros(solver.n_ks,solver.is_nt);
	else
		solver.snkt=zeros(solver.n_ks,solver.is_nfft);
		solver.sEkt=zeros(solver.n_ks,solver.is_nfft);
		solver.snkw=zeros(solver.n_ks,solver.is_nfft);
		solver.sEkw=zeros(solver.n_ks,solver.is_nfft);
		solver.is_fs=1/(solver.h*solver.is_stride)* ...
			[-solver.is_nfft/2:solver.is_nfft/2-1]/solver.is_nfft;
	end
	%solver.nks2=zeros(solver.n_ks,solver.is_nt);
	%solver.Eks2=zeros(solver.n_ks,solver.is_nt);
end
solver.n_is = 0;

