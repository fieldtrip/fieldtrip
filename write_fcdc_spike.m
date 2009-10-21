function write_fcdc_spike(varargin)

% this function is deprecated, please use the write_spike function instead

% Copyright (C) 2007-2009, Robert Oostenveld
%
% $Log: write_fcdc_spike.m,v $
% Revision 1.6  2009/05/07 14:21:29  roboos
% deprecated the read_fcdc and write_fcdc functions, give warning and mention the correct function to be used
%

fieldtripdefs

warning('this function is deprecated, please use the write_spike function instead');

% use the low-level writing function
write_spike(varargin{:});

