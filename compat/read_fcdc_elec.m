function [sens] = read_fcdc_elec(filename)

% this function is deprecated, please use the read_sens function instead

% Copyright (C) 2005-200, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

warning('this function is deprecated, please use the read_sens function instead');

% use the low-level reading function
[sens] = read_sens(filename);

