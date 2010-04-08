function [x, state] = preproc_standardize(x, begsample, endsample, state)

% PREPROC_STANDARDIZE performs a z-transformation or standardization
% of the data. The standardized data will have a zero-mean and a unit
% standard deviation.
%
% Use as
%   [dat] = preproc_standardize(dat, begsample, endsample)
% where
%   dat        data matrix (Nchans X Ntime)
%   begsample  index of the begin sample for the mean and stdev estimate
%   endsample  index of the end sample for the mean and stdev estimate
%
% If no begin and end sample are specified, it will be estimated on the
% complete data.
%
% See also PREPROC

% Copyright (C) 2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if nargin<2 || isempty(begsample)
  begsample = 1;
end

if nargin<3 || isempty(endsample)
  endsample = size(x,2);
end

if nargin<4
  state = [];
end

% get the data selection
y = x(:,begsample:endsample);

% determine the size of the selected data: nChans X nSamples
[m, n] = size(y);

% compute the sum and sum of squares
s  = sum(y,2);
ss = sum(y.^2,2);

% include the state information from the previous calls
if ~isempty(state)
  s  = s  + state.s;
  ss = ss + state.ss;
  n  = n  + state.n;
end

% compute the mean and standard deviation
my = s ./ n;
sy = sqrt((ss - (s.^2)./n) ./ (n-1));

% standardize the complete input data
x  = (x - repmat(my, 1, size(x, 2))) ./ repmat(sy, 1, size(x, 2));

% remember the state
state.s  = s;  % sum
state.ss = ss; % sum of sqares
state.n  = n;  % number of samples
