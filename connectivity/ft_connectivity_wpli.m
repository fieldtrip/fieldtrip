function [wpli, v, n] = ft_connectivity_wpli(input, varargin)

% FT_CONNECTIVITY_WPLI computes the weighted phase lag index from a
% data-matrix containing a cross-spectral density. It implements the method
% described in: Vinck M, Oostenveld R, van Wingerden M, Battaglia F,
% Pennartz CM. An improved index of phase-synchronization for
% electrophysiological data in the presence of volume-conduction, noise and
% sample-size bias. Neuroimage. 2011 Apr 15;55(4):1548-65.
%
% Use as
%   [wpi, v, n] = ft_connectivity_wpli(input, varargin)
% 
% The input data input should be organized as:
%
%   Repetitions x Channel x Channel (x Frequency) (x Time)
%
% or
%
%   Repetitions x Channelcombination (x Frequency) (x Time)
% 
% The first dimension should contain repetitions and should not contain an
% average already. Also, it should not consist of leave one out averages.
%
% Additional input arguments come as key-value pairs:
%
%   dojack   1 or 0,   compute a variance estimate, based on leave-one-out
%   feedback 'none', 'text', 'textbar' type of feedback showing progress of
%                   computation
%   debias 1 (or true) or 0 (or false), compute debiased wpli or not
%
% The output wpli contains the wpli, v is a leave-one-out variance estimate
% which is only computed if dojack = 1,and n is the number of repetitions
% in the input data.
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
debias      = ft_getopt(varargin, 'debias');
dojack      = ft_getopt(varargin, 'dojack');

siz = size(input);
n = siz(1);
ft_progress('init', feedback, 'computing metric...');
if n>1
  input    = imag(input);          % make everything imaginary  
  outsum   = nansum(input,1);      % compute the sum; this is 1 x size(2:end)
  outsumW  = nansum(abs(input),1); % normalization of the WPLI
  if debias
    outssq   = nansum(input.^2,1);
    wpli     = (outsum.^2 - outssq)./(outsumW.^2 - outssq); % do the pairwise thing in a handy way
  else
    wpli     = outsum./outsumW; % estimator of E(Im(X))/E(|Im(X)|)
  end    
  wpli = reshape(wpli,siz(2:end)); % remove the first singular dimension
else
  wpli = NaN(siz(2:end)); % for one observation, we should return NaNs
  ft_warning('ft_connectivity_wpli:nTrials', 'computation wpli requires >1 trial, returning NaNs');
end

[leave1outsum, leave1outssq] = deal(zeros([1 siz(2:end)]));
if dojack && n>2 % n needs to be larger than 2 to get a meaningful variance
  for k = 1:n
    s  = outsum  - input(k,:,:,:,:,:,:); % works for any array up to 7-D
    sw = outsumW - abs(input(k,:,:,:,:,:,:));
    if debias
      sq    = outssq - input(k,:,:,:,:,:,:).^2; 
      num   = s.^2  - sq;
      denom = sw.^2 - sq;
    else
      num   = s; % this is estimator of E(Im(X))
      denom = sw; % estimator of E(|Im(X)|)
    end        
    tmp          = num./denom; % avoids doing the division twice
    tmp(isnan(tmp)) = 0; % added for nan support
    leave1outsum = leave1outsum + tmp;% added this for nan support
    leave1outssq = leave1outssq + tmp.^2; % added this for nan support                              
  end  
  % compute the sem here 
  n = sum(~isnan(input),1); % this is the actual df when nans are found in the input matrix
  v = (n-1).^2.*(leave1outssq - (leave1outsum.^2)./n)./(n - 1); % 11.5 efron, sqrt and 1/n done in ft_connectivityanalysis
  v = reshape(v,siz(2:end)); % remove the first singular dimension   
  n = reshape(n,siz(2:end));  
elseif dojack && n<=2
  v = NaN(siz(2:end));
else
  v = [];
end
ft_progress('close');
