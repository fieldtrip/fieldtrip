function [spike] = read_fcdc_spike(filename)

% this function is deprecated, please use the read_spike function instead

% Copyright (C) 2007-2009, Robert Oostenveld
%
% $Log: read_fcdc_spike.m,v $
% Revision 1.5  2009/05/07 14:21:29  roboos
% deprecated the read_fcdc and write_fcdc functions, give warning and mention the correct function to be used
%

fieldtripdefs

warning('this function is deprecated, please use the read_spike function instead');

% use the low-level reading function
[spike] = read_spike(filename);

