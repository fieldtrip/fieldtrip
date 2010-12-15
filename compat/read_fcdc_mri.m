function [mri] = read_fcdc_mri(filename)

% this function is deprecated, please use the read_mri function instead

% Copyright (C) 2004-2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

ft_defaults

warning('this function is deprecated, please use the read_mri function instead');

% use the low-level reading function
[mri] = read_mri(filename);

