function [x, ut] = svdfft(f, n, trltapcnt);

% SVDFFT computes a rotated FFT matrix, using the real part of the
% cross-spectral density matrix. This rotation ensures that the phase
% relationship of the underlying sources does not change, while rotating
% the channels such that the first channel contains the maximal amplitude
% signal.
%
% Use as
%   [x, ut] = svdfft(f, n, trltapcnt);
% where
%   n           number of components (orientations) to keep in the output (e.g. 1 or 3)
%   trltapcnt   vector of length Ntrials with the number of tapers

% Copyright (C) 2005-2007, Robert Oostenveld & Jan-Mathijs Schoffelen
%
% $Log: svdfft.m,v $
% Revision 1.8  2008/04/21 14:32:38  jansch
% added the option to output a variable number of components, explaining n%
% of the variance.
%
% Revision 1.7  2007/01/17 17:04:43  roboos
% added a space, changed year
%
% Revision 1.6  2007/01/04 15:27:38  roboos
% updated documentation
%
% Revision 1.5  2006/05/08 09:00:19  roboos
% also return the rotation that is applied to the data
%
% Revision 1.4  2006/04/11 13:05:08  jansch
% fixed bug in loopvariable
%
% Revision 1.3  2006/03/24 14:56:41  jansch
% fixed typo
%
% Revision 1.2  2006/03/24 14:37:22  jansch
% included the option to equally weigh the trials in the case of unequal trial
% lengths. this needs an optional input argument trltapcnt. note that the output
% still contains the single taper estimates
%
% Revision 1.1  2005/09/09 08:20:47  roboos
% new implementation
%

if nargin == 1,
  n         = size(f,1);
  trltapcnt = ones(size(f,1),1);
elseif nargin == 2,
  trltapcnt = ones(size(f,1),1);
elseif nargin == 3,
  if isempty(n), n = size(f,1); end 
end

if all(trltapcnt==trltapcnt(1)),
  c = f * f';  
else
  trltapcnt = trltapcnt(:);
  sumtapcnt = cumsum([0;trltapcnt]);
  c         = zeros(size(f,1), size(f,1));
  for j = 1:length(sumtapcnt)-1
    c = c + [f(:, sumtapcnt(j)+1:sumtapcnt(j+1)) * f(:, sumtapcnt(j)+1:sumtapcnt(j+1))']./trltapcnt(j);
  end
end

if n==size(f,1),
  % do a complete decomposition
  [u, s, v] = svd(real(c));
elseif n<1,
  % do a complete decomposition and only keep the biggest components which together explain n percent of the variance 
  [u, s, v] = svd(real(c)); 
  s         = cumsum(diag(s))./sum(diag(s));
  n         = length(find(s<=n));
  u         = u(:, 1:n);
else
  % only decompose the first n components
  [u, s, v] = svds(real(c),n);
end

ut = u';      % this rotates the data in the direction of the maximum power
x  = ut * f;  % apply the rotation on the data

