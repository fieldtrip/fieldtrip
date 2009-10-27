function [event] = read_fcdc_event(filename)

% this function is deprecated, please use the read_event function instead

% Copyright (C) 2004-2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

warning('this function is deprecated, please use the read_event function instead');

% use the low-level reading function
[event] = read_event(filename);

