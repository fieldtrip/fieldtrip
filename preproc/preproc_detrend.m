function [dat, beta, x] = preproc_detrend(dat, begsample, endsample, order)

% PREPROC_DETREND removes linear or higher order polynomial trends from the
% data using using General Linear Modeling
%
% Use as
%   [dat] = preproc_detrend(dat, begin, end, order)
% where
%   dat        data matrix (Nchans X Ntime)
%   begsample  index of the begin sample for the trend estimate
%   endsample  index of the end sample for the trend estimate
%   order      number representing the polynomial order (default = 1, i.e. linear)
%
% If no begin and end sample are specified for the trend estimate, it
% will be estimated on the complete data.
%
% See also PREPROC

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: preproc_detrend.m,v $
% Revision 1.3  2008/06/10 16:03:40  roboos
% fixed small bug and typ, thanks to Saskia
%
% Revision 1.2  2008/05/23 09:13:58  roboos
% cleaned up code and documentation, ensure that all functions are consistent, added proper implementation to the scratch functions
%

% determine the size of the data
[Nchans, Nsamples] = size(dat);

% determine the interval to use for baseline correction
if nargin<2 || isempty(begsample)
  begsample = 1;
end
if nargin<3 || isempty(endsample)
  endsample = Nsamples;
end

% determine the order of the polynomial trend to be removed, default is linear
if nargin<4 || isempty(order)
  order = 1;
end

% create a matrix with regressor components
basis    = 1:Nsamples;
x        = zeros(order+1,Nsamples);
for i=0:order
  x(i+1,:) = basis.^(i);
end
% estimate the polynomial trend using General Linear Modeling, where dat=beta*x+noise
beta = dat(:,begsample:endsample)/x(:,begsample:endsample);
% subtract the trend from the complete data
dat  = dat - beta*x;
