function [hdr] = read_fcdc_header(filename)

% this function is deprecated, please use the read_header function instead

% Copyright (C) 2003-2009, Robert Oostenveld
%
% $Log: read_fcdc_header.m,v $
% Revision 1.49  2009/05/07 14:21:29  roboos
% deprecated the read_fcdc and write_fcdc functions, give warning and mention the correct function to be used
%

fieldtripdefs

warning('this function is deprecated, please use the read_header function instead');

% use the low-level reading function
[hdr] = read_header(filename);

