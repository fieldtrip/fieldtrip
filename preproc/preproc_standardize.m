function x = preproc_standardize(x, begsample, endsample)

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
% $Log: preproc_standardize.m,v $
% Revision 1.2  2009/09/30 13:02:03  jansch
% fixed typo in comment
%
% Revision 1.1  2008/11/10 20:58:32  roboos
% new function, original was used in trialfun_emg
%

% determine the size of the input data: nChans X nSamples
[m,n] = size(x);

% this function operates along the 2nd dimension
dim   = 2;

if nargin>1
  mx = mean(x(:,begsample:endsample),dim);
  sx = std(x(:,begsample:endsample),0,dim);
else
  mx = mean(x,dim);
  sx = std(x,0,dim);
end

x     = (x - repmat(mx,[1 n]))./repmat(sx,[1 n]);

