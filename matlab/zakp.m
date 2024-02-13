function s=param()
% function s=param()

%
% $Id: zakp.m,v 1.9 2011/03/26 09:20:42 patrick Exp $
%
% Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

s.n          = 4096;          % number of components
s.mmax       = 1365;          % samples kept to zero after mode mmax
s.n2         = fix(s.n/2);
s.nr         = 126;           % sample no seen by the radar
s.mass_ratio = 40104;         % M/m
s.t_ratio    = 2.0;           % Te/Ti
s.ecoll      = 1.4;           % electron collision freq
s.nyic       = 0;             % ion collisional damping, F region
s.domega     = 25;            % frequency mismatch
s.h          = 1e-4;          % time-step length
s.h2         = s.h/2;
s.im         = sqrt(-1); 
s.length     = 2*pi/0.0984;   % length of simulation cell
s.nstart     = 200000;        % starting timestep for ave.
s.l2         = s.length/2;
s.dx         = s.length/s.n;
s.pi2        = 2*pi;
s.ndx        = s.n*s.dx;
s.pi2ndx     = s.pi2/s.ndx;   % distance between k-components
s.pi2ndx2    = s.pi2ndx^2;
s.edrive     = 1.0;           % electric pump field
s.fact1      = 2e-3;          % initial E field fluctuation
s.fact2      = 0e-3;          % initial ion fluctuation
s.profile    = 2;             % 1 is envelope soliton, 2 is random fluctuations
s.k2s        = s.pi2ndx2*(s.n-s.mmax+1)^2; % mode m = n-mmax+1
s.ftn        = 16384;         % number of samples for freq spectra
s.nbreak     = 1020;          % breaking point for landau damping
s.choice     = 1;
s.nq         = 102400;
s.seed       = 1;
s.t_out      = 102400;        % flag for |E|^2 and n plots

% soliton params
s.v          = .6;            % soliton velocity
s.omega      = 0;             % carrier frequency
s.x0         = 20;            % distance from x=0

% spectra
s.ns=160;

% ode solver
s.k = ifftshift(2*pi*[-s.n2:s.n2-1]'/s.length);
%s.k(1)
s.ik = i*s.k;
s.k2 = s.k.^2;
s.ik2 = i*s.k2;

% zero padded aliased modes
s.zm = [s.mmax+2:s.n+1-s.mmax-1]';
% computed modes [k=0, k>0, k<0]
s.cm = [1, 1+[1:s.mmax], s.n+1-fliplr([1:s.mmax])]';

% number of computed modes
s.ncm = length(s.cm);
% xlim for zero-centered spectra 
s.klim = s.n/2+1+[-s.mmax-1, s.mmax+1]';
if s.ncm+2 ~= s.klim(2)-s.klim(1)+1,
  s.ncm
  s.klim(2)-s.klim(1)+1
  error('klim wrongly defined');
end

s.odesteps = fix(s.nq/s.ns);
s.t0 = 0;


