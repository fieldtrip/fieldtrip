function write_fcdc_spike(varargin)

% this function is deprecated, please use the write_spike function instead

% Copyright (C) 2007-2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

warning('this function is deprecated, please use the write_spike function instead');

% use the low-level writing function
write_spike(varargin{:});

