function [dtf, dtfvar, n] = ft_connectivity_dtf(input, varargin)

% FT_CONNECTIVITY_DTF computes directed transfer function.
% 
% Use as
%   [d, v, n] = ft_connectivity_dtf(h, key1, value1, ...)
%
% Input arguments: 
%   h = spectral transfer matrix, Nrpt x Nchan x Nchan x Nfreq (x Ntime),
%      Nrpt can be 1.
%
% additional options need to be specified as key-value pairs and are:
%   'hasjack'  = 0 (default) is a boolean specifying whether the input
%                contains leave-one-outs, required for correct variance
%                estimate.
%   'feedback' = string, determining verbosity (default = 'none'), see FT_PROGRESS
%   'crsspctrm' = matrix containing the cross-spectral density. If this 
%                 matrix is defined, the function
%                 returns the ddtf, which requires an estimation of partial
%                 coherence from this matrix.
%   'invfun'   = 'inv' (default) or 'pinv', the function used to invert the
%                crsspctrm matrix to obtain the partial coherence. Pinv is
%                useful if the data are poorly-conditioned.
%
%
% Output arguments:
%   d = partial directed coherence matrix Nchan x Nchan x Nfreq (x Ntime).
%       If multiple observations in the input, the average is returned.
%   v = variance of d across observations.
%   n = number of observations.
%
% Typically, nrpt should be 1 (where the spectral transfer matrix is
% computed across observations. When nrpt>1 and hasjack is true the input
% is assumed to contain the leave-one-out estimates of H, thus a more
% reliable estimate of the relevant quantities.

% Copyright (C) 2009-2017, Jan-Mathijs Schoffelen
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

hasjack   = ft_getopt(varargin, 'hasjack', 0);
powindx   = ft_getopt(varargin, 'powindx');
feedback  = ft_getopt(varargin, 'feedback', 'none');
crsspctrm = ft_getopt(varargin, 'crsspctrm');
invfun    = ft_getopt(varargin, 'invfun', 'inv');

switch invfun
  case {'inv' 'pinv'}
    invfun = str2func(invfun);
  otherwise
    ft_error('unknown specification of inversion-function for the transfer matrix');
end

if ~isempty(powindx)
  % this error message is rather uninformative, but is kept for now for
  % backward compatibility reasons (i.e. it might exist when called from
  % ft_connectivityanalysis
  ft_error('linearly indexed data for dtf computation is at the moment not supported');
end

siz    = [size(input) 1];
n      = siz(1);
outsum = zeros(siz(2:end));
outssq = zeros(siz(2:end));

% check the crsspctrm, if it's present, compute the partial coherence
if ~isempty(crsspctrm)
  assert(isequal(size(crsspctrm),size(input)), 'the input data should be of the same size as the crsspctrm');
  fprintf('computing dDTF in the presence of a crsspctrm\n');
  
  % the crsspctrm allows for the partial coherence to be computed
  pdim   = prod(siz(4:end));
  tmpcrs = reshape(crsspctrm, [siz(1:3) pdim]);
  ft_progress('init', feedback, 'computing partial coherence...');
  for k = 1:n
    ft_progress(k/n, 'computing partial coherence for replicate %d from %d\n', k, n);
    tmp = reshape(tmpcrs(k,:,:,:), [siz(2:3) pdim]);
    for m = 1:pdim
      tmp(:,:,m) = invfun(tmp(:,:,m));
      tmp(:,:,m) = abs(tmp(:,:,m))./sqrt(abs(diag(tmp(:,:,m))*diag(tmp(:,:,m))'));
    end
    tmpcrs(k,:,:,:) = tmp;
  end
  ft_progress('close');
  crsspctrm = reshape(tmpcrs, siz);
end

% data should be represented as chan_chan_therest
for j = 1:n
  tmph   = reshape(input(j,:,:,:,:), siz(2:end));
  if isempty(crsspctrm)
    % plain DTF
    den    = sum(abs(tmph).^2,2);
    tmpdtf = abs(tmph)./sqrt(repmat(den, [1 siz(2) 1 1 1]));
  else
    % dDTF
    tmpc   = reshape(crsspctrm(j,:,:,:,:), siz(2:end));
    
    den    = sum(sum(abs(tmph).^2,3),2);
    tmpdtf = abs(tmph)./sqrt(repmat(den, [1 siz(2) siz(4) 1 1 1]));
    tmpdtf = tmpdtf.*tmpc;
  end
  outsum = outsum + tmpdtf;
  outssq = outssq + tmpdtf.^2;
end
dtf = outsum./n;

if n>1, %FIXME this is strictly only true for jackknife, otherwise other bias is needed
  if hasjack,
    bias = (n - 1).^2;
  else
    bias = 1;
  end
  dtfvar = bias.*(outssq - (outsum.^2)/n)./(n-1);
else
  dtfvar = [];
end

% this is to achieve the same convention for all connectivity metrics:
% row -> column
for k = 1:prod(siz(4:end))
  dtf(:,:,k)    = transpose(dtf(:,:,k));
  if ~isempty(dtfvar)
    dtfvar(:,:,k) = transpose(dtfvar(:,:,k));
  end
end
