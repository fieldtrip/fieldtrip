function [c] = ft_connectivity_powcorr_ortho(mom, varargin)

% FT_CONNECTIVITY_POWCORR_ORTHO computes power correlation after removing
% the zero-lag contribution on a trial-by-trial basis, according to Hipp's
% Nature Neuroscience paper. 
%
% Use as
%   c = ft_connectivity_powcorr(mom)
%   c = ft_connectivity_powcorr(mom, 'refindx', refindx)
%
% Where mom is a NchanxNrpt matrix containing the complex-valued amplitude
% and phase information at a given frequency, and the optional key refindx
% specifies the index/indices of the channels that serve as a reference 
% channel. (Default is 'all').
%
% The output c is a NchanxNrefchan matrix that contain the power correlation
% for all channels orthogonalised relative to the reference channel in the first
% Nrefchan columns, and the power correlation for the reference channels 
% orthogonalised relative to the channels in the second Nrefchan columns.

% Copyright (C) 2012 Jan-Mathijs Schoffelen
% 2021: Andrea Ibarra-Chaoul + Tobias Ludwig:
% added feature to handle multiple tapers by computing orth. power-correlation
% for each taper separately and then averaging over tapers
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$


refindx = ft_getopt(varargin, 'refindx', 'all');
tapvec  = ft_getopt(varargin, 'tapvec',  ones(1,size(mom,2)));
% tapvec = [ntap, ntap, ...], same number of tapers for each rep. = trial

if strcmp(refindx, 'all')
  refindx = 1:size(mom,1);
end

nchan = size(mom,1);
ntap  = tapvec(1);
nrpt  = numel(tapvec); % number of trials / repetitions

if ~all(tapvec==ntap)
  ft_error('unequal number of tapers per observation is not yet supported');
end

% final connectivity matrix
c = zeros(nchan, numel(refindx))+nan;

% only need to do these two things once (out of next forloop)
cXnorm = conj(mom./abs(mom));
powX   = abs(mom).^2;

for k = 1:numel(refindx) % for each source/channel
  indx   = refindx(k);
  target = setdiff(1:size(mom,1), indx);
  
  Y    = repmat(mom(indx,:), [nchan, 1]); % Y = y nchan times stacked
   
  %% orthogonalization in one direction: Y wrt X
  powYorth = abs(imag(Y.*cXnorm)).^2;
  
  zYorth   = zeros(nchan, nrpt*ntap);
  zX       = zeros(nchan, nrpt*ntap);
  
  % correlation for each taper separately
  for tap = 1:ntap
    idx = tap + ntap * (0:(nrpt-1));
    zYorth(target,idx) = standardise(log10(powYorth(target,idx)), 2);
    zX(target,idx)     = standardise(log10(powX(target,idx)), 2);
  end
  
  c1 = mean(zX.*zYorth, 2); % take correlation averaging over trials+tapers
  
  %% in the other direction: orthogonalize X wrt Y
  cYnorm = conj(Y./abs(Y));
  
  powXorth = abs(imag(mom.*cYnorm)).^2;
  powY     = abs(Y).^2;
  
  zXorth   = zeros(nchan, nrpt*ntap);
  zY       = zeros(nchan, nrpt*ntap);
  
  for tap = 1:ntap
    idx = tap + ntap * (0:(nrpt-1));
    zXorth(target,idx) = standardise(log10(powXorth(target,idx)), 2);
    zY(target,idx)     = standardise(log10(powY(target,idx)), 2);
  end
  
  c2 = mean(zXorth.*zY, 2);

  c(:,k) = (c1+c2)./2;
  
end
end
