function [pdc, pdcvar, n] = ft_connectivity_pdc(input, varargin)

% FT_CONNECTIVITY_PDC computes partial directed coherence.
% 
% Use as
%   [p, v, n] = ft_connectivity_pdc(h, key1, value1, ...)
%
% Input arguments: 
%   H = spectral transfer matrix, Nrpt x Nchan x Nchan x Nfreq (x Ntime),
%      Nrpt can be 1.
%
% additional options need to be specified as key-value pairs and are:
%   'hasjack'  = 0 (default) is a boolean specifying whether the input
%                contains leave-one-outs, required for correct variance
%                estimate
%   'feedback' = string, determining verbosity (default = 'none'), see FT_PROGRESS
%   'invfun'   = 'inv' (default) or 'pinv', the function used to invert the
%                transfer matrix to obtain the fourier transform of the
%                MVAR coefficients. Use 'pinv' if the data are
%                poorly-conditioned.
%   'noisecov' = matrix containing the covariance of the residuals of the
%                MVAR model. If this matrix is defined, the function
%                returns the generalized partial directed coherence.
%
% Output arguments:
%   p = partial directed coherence matrix Nchan x Nchan x Nfreq (x Ntime).
%       If multiple observations in the input, the average is returned.
%   v = variance of p across observations.
%   n = number of observations.
%
% Typically, nrpt should be 1 (where the spectral transfer matrix is
% computed across observations. When nrpt>1 and hasjack is true the input
% is assumed to contain the leave-one-out estimates of H, thus a more
% reliable estimate of the relevant quantities.
%
% This function implements the metrices described in:
%  - Baccala et al., Biological Cybernetics 2001, 84(6), 463-74.
%  - Baccala et al., 15th Int.Conf.on DSP 2007, 163-66.
%
% The implemented algorithm has been tested against the implementation in
% the SIFT-toolbox. It yields numerically identical results to what is
% known in the SIFT-toolbox as 'nPDC' (for PDC) and 'GPDC' for generalized
% pdc.

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

hasjack  = ft_getopt(varargin, 'hasjack', 0);
feedback = ft_getopt(varargin, 'feedback', 'none');
invfun   = ft_getopt(varargin, 'invfun',   'inv');
noisecov = ft_getopt(varargin, 'noisecov');

switch invfun
  case {'inv' 'pinv'}
    invfun = str2func(invfun);
  otherwise
    ft_error('unknown specification of inversion-function for the transfer matrix');
end

% crossterms are described by chan_chan_therest
siz = [size(input) 1];
n   = siz(1);

if ~isempty(noisecov)
  ft_progress('init', feedback, 'computing generalized partial directed coherence...');
  scale = diag(sqrt(1./diag(noisecov))); % get the 1./sqrt(var) of the signals
else
  ft_progress('init', feedback, 'computing partial directed coherence...');
  scale = eye(siz(2));
end

% pre-allocate some variables
outsum = zeros(siz(2:end));
outssq = zeros(siz(2:end));

% the mathematics for pdc is most straightforward using the inverse of the
% transfer function
pdim     = prod(siz(4:end));
tmpinput = reshape(input, [siz(1:3) pdim]);
ft_progress('init', feedback, 'inverting the transfer function...');
for k = 1:n
  ft_progress(k/n, 'inverting the transfer function for replicate %d from %d\n', k, n);
  tmp = reshape(tmpinput(k,:,:,:), [siz(2:3) pdim]);
  for m = 1:pdim
    tmp(:,:,m) = scale*invfun(tmp(:,:,m));
  end
  tmpinput(k,:,:,:) = tmp;
end
ft_progress('close');
input = reshape(tmpinput, siz);

for j = 1:n
  ft_progress(j/n, 'computing metric for replicate %d from %d\n', j, n);
  invh   = reshape(input(j,:,:,:,:), siz(2:end));
  
  
  den    = sum(abs(invh).^2,1);
  tmppdc = abs(invh)./sqrt(repmat(den, [siz(2) 1 1 1 1]));
  
  
  outsum = outsum + tmppdc;
  outssq = outssq + tmppdc.^2;
end
ft_progress('close');

pdc = outsum./n;

if n>1,
  if hasjack
    bias = (n-1).^2;
  else
    bias = 1;
  end
  pdcvar = bias*(outssq - (outsum.^2)./n)./(n - 1);
else
  pdcvar = [];
end

% this is to achieve the same convention for all connectivity metrics:
% row -> column
for k = 1:prod(siz(4:end))
  pdc(:,:,k)    = transpose(pdc(:,:,k));
  if ~isempty(pdcvar)
    pdcvar(:,:,k) = transpose(pdcvar(:,:,k));
  end
end
