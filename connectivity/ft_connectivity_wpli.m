function [wpli, v, n] = ft_connectivity_wpli(input, varargin)

% FT_CONNECTIVITY_WPLI computes WPLI from a data-matrix
% containing a cross-spectral density
%
% Use as
%   [wpi, v, n] = ft_connectivity_wpli(input, varargin)
% 
% The input data input should be organized as:
%   Repetitions x Channel x Channel (x Frequency) (x Time)
% or
%   Repetitions x Channelcombination (x Frequency) (x Time)
% 
% The first dimension should contain repetitions and should not contain an average already.
% Also, it should not consist of leave one out averages.
%
% Additional input arguments come as key-value pairs:
%
% feedback 'none', 'text', 'textbar' type of feedback showing progress of
%                   computation
% debias 1 (or true) or 0 (or false), we compute wpli or debiased wpli (Vinck et al., 2011)
%
% The output wpli contains the wpli, v is a leave-one-out variance estimate
% which is only computed if dojack = 1,and n is the number of repetitions in the input data.
% 
% This is a helper function to FT_CONNECTIVITYANALYSIS
% 
% See also FT_CONNECTIVITYANALYSIS

% Copyright (C) 2011, Martin Vinck 
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

feedback    = keyval('feedback', varargin{:}); if isempty(feedback), feedback = 'none'; end
debias      = keyval('debias',   varargin{:});
dojack      = keyval('dojack',   varargin{:});

siz = size(input);
n = siz(1);
ft_progress('init', feedback, 'computing metric...');
if n>1
  input    = imag(input);        % make everything imaginary  
  outsum   = nansum(input);      % compute the sum; this is 1 x size(2:end)
  outsumW  = nansum(abs(input)); % normalization of the WPLI
  if debias
    outssq   = nansum(input.^2);
    wpli     = (outsum.^2 - outssq)./(outsumW.^2 - outssq); % do the pairwise thing in a handy way
  else
    wpli     = outsum./outsumW; % estimator of E(Im(X))/E(|Im(X)|)
  end    
  wpli = reshape(wpli,siz(2:end)); % remove the first singular dimension
else
  wpli = NaN(siz(2:end)); % for one observation, we should return NaNs
  warning('ft_connectivity_wpli:nTrials', 'computation wpli requires >1 trial, returning NaNs');
end

[leave1outsum, leave1outssq] = deal(0);
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
