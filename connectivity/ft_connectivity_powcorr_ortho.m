function [c] = ft_connectivity_powcorr_ortho(inputdata, varargin)

% FT_CONNECTIVITY_POWCORR_ORTHO computes power correlation after removing
% the zero-lag contribution on a trial-by-trial basis, according to Hipp's
% Nature Neuroscience paper.
%
% Use as
%   [c] = ft_connectivity_powcorr(inputdata, ...)
%
% Where the input is a Nchan*Nrpt matrix containing the complex-valued amplitude
% and phase information at a given frequency.
%
% The output c is a Nchan*Nref matrix that contain the power correlation for all
% channels orthogonalised relative to the reference channel in the first Nref
% columns, and the power correlation for the reference channels orthogonalised
% relative to the channels in the second Nref columns.
%
% Additional optional input arguments come as key-value pairs:
%   'refindx'  = index/indices of the channels that serve as a reference channel (default is all)
%   'tapvec'   = vector with the number of tapers per trial
%
% See also CONNECTIVITY, FT_CONNECTIVITYANALYSIS

% Copyright (C) 2012 Jan-Mathijs Schoffelen
% Copyright (C) 2021 Andrea Ibarra-Chaoul and Tobias Ludwig, who added the feature to
% handle multiple tapers by computing orthogonal power-correlation for each taper
% separately and then averaging over tapers
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
tapvec  = ft_getopt(varargin, 'tapvec',  ones(1,size(inputdata,2))); % default is 1 taper per trial

if strcmp(refindx, 'all')
  refindx = 1:size(inputdata,1);
end

[nchan, nrpttap] = size(inputdata);
ntap  = tapvec(1);
nrpt  = numel(tapvec); % number of trials / repetitions

if sum(tapvec)~=nrpttap
  ft_error('the number of tapers and trials does not match');
end

if ~all(tapvec==ntap)
  ft_error('unequal number of tapers per observation is not yet supported');
end

% final connectivity matrix
c = zeros(nchan, numel(refindx))+nan;

% only need to do these two things once (out of next forloop)
cXnorm = conj(inputdata./abs(inputdata));
powX   = abs(inputdata).^2;

for k = 1:numel(refindx) % for each source/channel
  indx   = refindx(k);
  target = setdiff(1:size(inputdata,1), indx);
  
  Y    = repmat(inputdata(indx,:), [nchan, 1]); % Y = y nchan times stacked
  
  % orthogonalization in one direction: Y wrt X
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
  
  % in the other direction: orthogonalize X wrt Y
  cYnorm = conj(Y./abs(Y));
  
  powXorth = abs(imag(inputdata.*cYnorm)).^2;
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
