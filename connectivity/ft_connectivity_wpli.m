function [wpli, v, n] = ft_connectivity_wpli(inputdata, varargin)

% FT_CONNECTIVITY_WPLI computes the weighted phase lag index from a data matrix
% containing the cross-spectral density. This implements the method described in
% Vinck M, Oostenveld R, van Wingerden M, Battaglia F, Pennartz CM. An improved index
% of phase-synchronization for electrophysiological data in the presence of
% volume-conduction, noise and sample-size bias. Neuroimage. 2011 Apr
% 15;55(4):1548-65.
%
% Use as
%   [wpi, v, n] = ft_connectivity_wpli(inputdata, ...)
%
% The input data input should contain cross-spectral densities organized as:
%   Repetitions x Channel x Channel (x Frequency) (x Time)
% or
%   Repetitions x Channelcombination (x Frequency) (x Time)
% 
% Alternatively, the input data can contain fourier coefficients organized
% as:
%   Repetitions_tapers x Channel (x Frequency) (x Time) 
%
% The first dimension of the input data matrix should contain repetitions and should not 
% contain an average already. Also, the input should not consist of leave-one-out averages.
%
% The output wpli contains the wpli, v is a leave-one-out variance estimate
% which is only computed if dojack=true, and n is the number of repetitions
% in the input data.
%
% Additional optional input arguments come as key-value pairs:
%   'dojack'    = boolean, compute a variance estimate based on
%                 leave-one-out, only supported when input data is a
%                 bivariate cross-spectral density
%   'debias'    = boolean, compute debiased wpli or not
%   'feedback'  = 'none', 'text', 'textbar', 'dial', 'etf', 'gui' type of feedback 
%                 showing progress of computation, see FT_PROGRESS
%   'isunivariate' = boolean, compute CSD on fly (saves memory with many trials)
%   'cumtapcnt' = vector that contains the cumulative taper counter, defining how
%                 tapers should be combined to define repetitions. If not
%                 defined (or empty), it will be ones(size(input,1),1),
%                 i.e. each slice of the matrix is considered a repetition.
%                 This option is only function in case isunivariate = true
%
% See also FT_CONNECTIVITYANALYSIS

% Copyright (C) 2011-2022, Martin Vinck and Jan-Mathijs Schoffelen
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

feedback     = ft_getopt(varargin, 'feedback',     'none');
debias       = ft_getopt(varargin, 'debias');
dojack       = ft_getopt(varargin, 'dojack',       false);
isunivariate = ft_getopt(varargin, 'isunivariate', false);
cumtapcnt    = ft_getopt(varargin, 'cumtapcnt',    []);

if dojack && isunivariate
  error('jackknife variance estimates with on-the-fly csd computation is not supported');
end
if isunivariate
  if isempty(cumtapcnt)
    cumtapcnt = ones(size(inputdata,1), 1);
  end
  assert(sum(cumtapcnt)==size(inputdata,1));

  siz    = [size(inputdata) 1]; 
  nchan  = siz(2);
  outsiz = [nchan nchan siz(3:end)];
  n      = size(cumtapcnt,1);
  
  sumtapcnt = [0;cumsum(cumtapcnt(:,1))];

  outsum  = complex(zeros(outsiz));
  outsumW = complex(zeros(outsiz));
  outssq  = complex(zeros(outsiz));

  for k = 1:n 
    indx  = (sumtapcnt(k)+1):sumtapcnt(k+1);
    for m = 1:prod(outsiz(3:end))
      trial = transpose(inputdata(indx,:,m));
      csdimag = imag(trial*trial')./length(indx);

      outsum(:,:,m)  = outsum(:,:,m)  + csdimag;
      outsumW(:,:,m) = outsumW(:,:,m) + abs(csdimag);
      outssq(:,:,m)  = outssq(:,:,m)  + (csdimag.^2);
    end
  end
  if debias
    wpli = (outsum.^2 - outssq)./(outsumW.^2 - outssq);
  else
    wpli = outsum./outsumW;
  end
  v = [];
else
  siz = [ size(inputdata) 1 ];
  n = siz(1);
  if n>1
    inputdata    = imag(inputdata);          % make everything imaginary
    outsum   = nansum(inputdata,1);      % compute the sum; this is 1 x size(2:end)
    outsumW  = nansum(abs(inputdata),1); % normalization of the WPLI
    if debias
      outssq   = nansum(inputdata.^2,1);
      wpli     = (outsum.^2 - outssq)./(outsumW.^2 - outssq); % do the pairwise thing in a handy way
    else
      wpli     = outsum./outsumW; % estimator of E(Im(X))/E(|Im(X)|)
    end
    wpli = reshape(wpli,siz(2:end));
  else
    wpli = NaN(siz(2:end)); % for one observation, we should return NaNs
    ft_warning('ft_connectivity_wpli:nTrials', 'computation wpli requires >1 trial, returning NaNs');
  end

  [leave1outsum, leave1outssq] = deal(zeros([1 siz(2:end)]));
  if dojack && n>2 % n needs to be larger than 2 to get a meaningful variance
    ft_progress('init', feedback, 'computing metric...');
    for k = 1:n
      ft_progress(k/n, 'computing metric for replicate %d from %d\n', k, n);
      s  = outsum  - inputdata(k,:,:,:,:,:,:); % works for any array up to 7-D
      sw = outsumW - abs(inputdata(k,:,:,:,:,:,:));
      if debias
        sq    = outssq - inputdata(k,:,:,:,:,:,:).^2;
        num   = s.^2  - sq;
        denom = sw.^2 - sq;
      else
        num   = s; % this is estimator of E(Im(X))
        denom = sw; % estimator of E(|Im(X)|)
      end
      tmp          = num./denom;            % avoids doing the division twice
      tmp(isnan(tmp)) = 0;                  % added for nan support
      leave1outsum = leave1outsum + tmp;    % added this for nan support
      leave1outssq = leave1outssq + tmp.^2; % added this for nan support
    end
    ft_progress('close');

    % compute the sem here
    n = sum(~isnan(inputdata),1); % this is the actual df when nans are found in the input matrix
    v = (n-1).^2.*(leave1outssq - (leave1outsum.^2)./n)./(n - 1); % 11.5 efron, sqrt and 1/n done in ft_connectivityanalysis
    v = reshape(v,siz(2:end)); % remove the first singular dimension
    n = reshape(n,siz(2:end));
  elseif dojack && n<=2
    v = NaN(siz(2:end));
  else
    v = [];
  end
end
