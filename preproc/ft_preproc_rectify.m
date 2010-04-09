function [dat] = ft_preproc_rectify(dat)

% FT_PREPROC_RECTIFY rectifies the data, i.e. converts all samples with a
% negative value into the similar magnitude positive value
%
% Use as
%   [dat] = ft_preproc_rectify(dat)
% where
%   dat        data matrix (Nchans X Ntime)
%
% See also PREPROC

% Copyright (C) 2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

dat = abs(dat);
