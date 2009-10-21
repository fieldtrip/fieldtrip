function dat = preproc_medianfilter(dat, order);

% PREPROC_MEDIANFILTER applies a median filter, which smooths the data with
% a boxcar-like kernel except that it keeps steps in the data. This
% function requires the Matlab Signal Processing toolbox.
%
% Use as
%   [dat] = preproc_medianfilter(dat, order)
% where
%   dat        data matrix (Nchans X Ntime)
%   order      number, the length of the median filter kernel (default = 25)
%
% See also PREPROC

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: preproc_medianfilter.m,v $
% Revision 1.2  2008/05/23 09:13:58  roboos
% cleaned up code and documentation, ensure that all functions are consistent, added proper implementation to the scratch functions
%

% set the default filter order
if nargin<2 || isempty(order)
   order = 25;
end

dat = medfilt1(dat, order, [], 2);
