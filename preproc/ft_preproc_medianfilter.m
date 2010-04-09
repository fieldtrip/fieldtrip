function dat = ft_preproc_medianfilter(dat, order);

% FT_PREPROC_MEDIANFILTER applies a median filter, which smooths the data with
% a boxcar-like kernel except that it keeps steps in the data. This
% function requires the Matlab Signal Processing toolbox.
%
% Use as
%   [dat] = ft_preproc_medianfilter(dat, order)
% where
%   dat        data matrix (Nchans X Ntime)
%   order      number, the length of the median filter kernel (default = 25)
%
% See also PREPROC

% Copyright (C) 2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% set the default filter order
if nargin<2 || isempty(order)
   order = 25;
end

dat = medfilt1(dat, order, [], 2);
