function [c, v, n] = ft_connectivity_ppc(input, varargin)

% FT_CONNECTIVITY_PPC computes pairwise phase consistency or weighted pairwise phase
% consistency from a data-matrix containing a cross-spectral density. This implements
% the method described in Vinck M, van Wingerden M, Womelsdorf T, Fries P, Pennartz
% CM. The pairwise phase consistency: a bias-free measure of rhythmic neuronal
% synchronization. Neuroimage. 2010 May 15;51(1):112-22.
%
% Use as
%   [c, v, n] = ft_connectivity_ppc(input, ...)
%
% The input data input should be organized as:
%   Repetitions x Channel x Channel (x Frequency) (x Time)
% or
%   Repetitions x Channelcombination (x Frequency) (x Time)
%
% The first dimension should contain repetitions and should not contain an average
% already. Also, it should not consist of leave-one-out averages.
%
% Additional optional input arguments come as key-value pairs:
%   feedback  = 'none', 'text', 'textbar' type of feedback showing progress  of computation
%   weighted  = 1 (or true) or 0 (or false), we compute unweighted ppc or
%               weighted ppc, the weighting is according to the magnitude of
%               the cross-spectrum
%
% The output c contains the ppc, v is a leave-one-out variance estimate which is only
% computed if dojack = 1,and n is the number of repetitions in the input data.
%
% See also FT_CONNECTIVITYANALYSIS

% Copyright (C) 2011, Martin Vinck
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

feedback    = ft_getopt(varargin, 'feedback', 'none');
weighted    = ft_getopt(varargin, 'weighted');
dojack      = ft_getopt(varargin, 'dojack');

siz = size(input);
n = siz(1);
ft_progress('init', feedback, 'computing metric...');
if n>1
  if ~weighted
    input    = input./abs(input);  % normalize the crosspectrum
    outsum   = nansum(input);      % compute the sum; this is 1 x size(2:end)
    c        = (outsum.*conj(outsum) - n)./(n*(n-1)); % do the pairwise thing in a handy way
  else
    outsum   = nansum(input); % normalization of the WPLI
    outssq   = nansum(input.*conj(input));
    outsumw  = nansum(abs(input));
    c        = (outsum.*conj(outsum) - outssq)./(outsumw.*conj(outsumw) - outssq); % do the pairwise thing in a handy way
  end
  c          = reshape(c,siz(2:end)); % remove the first singular dimension
else
  c = NaN(siz(2:end)); % for one observation, we should return NaNs
  ft_warning('computation ppc requires >1 trial, returning NaNs')
end

[leave1outsum, leave1outssq] = deal(0);
if dojack && n>2 % n needs to be larger than 2 to get a meaningful variance
  for k = 1:n
    % this code works with both formats of input, also if it is 5-D
    s       = outsum - input(k,:,:,:,:,:,:); % index up to 7-D, this also works for 5-D then.
    if ~weighted
      num   = s.*conj(s) - (n-2);
      denom = (n-1)*(n-2);
    else
      sq    = outssq  - input(k,:,:,:,:,:,:).*conj(input(k,:,:,:,:,:,:));
      sw    = outsumw - abs(input(k,:,:,:,:,:,:));
      num   = s.*conj(s)   - sq;
      denom = sw.*conj(sw) - sq;
    end
    leave1outsum = leave1outsum + num./denom;
    leave1outssq = leave1outssq + (num./denom).^2;
  end
  % compute the sem here
  v = (n-1).^2*(leave1outssq - (leave1outsum.^2)./n)./(n - 1); % 11.5 efron, sqrt and 1/n done in ft_connectivityanalysis
  v = reshape(v,siz(2:end)); % remove the first singular dimension
elseif dojack && n<=2
  v = NaN(siz(2:end));
else
  v = [];
end
ft_progress('close');
