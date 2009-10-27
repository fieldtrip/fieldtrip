function [spike] = read_fcdc_spike(filename)

% this function is deprecated, please use the read_spike function instead

% Copyright (C) 2007-2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

warning('this function is deprecated, please use the read_spike function instead');

% use the low-level reading function
[spike] = read_spike(filename);

