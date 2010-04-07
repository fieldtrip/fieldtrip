function [hdr] = read_fcdc_header(filename)

% this function is deprecated, please use the read_header function instead

% Copyright (C) 2003-2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

warning('this function is deprecated, please use the read_header function instead');

% use the low-level reading function
[hdr] = read_header(filename);

