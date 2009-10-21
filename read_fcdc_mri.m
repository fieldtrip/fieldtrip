function [mri] = read_fcdc_mri(filename)

% this function is deprecated, please use the read_mri function instead

% Copyright (C) 2004-2009, Robert Oostenveld
%
% $Log: read_fcdc_mri.m,v $
% Revision 1.21  2009/05/07 14:21:29  roboos
% deprecated the read_fcdc and write_fcdc functions, give warning and mention the correct function to be used
%

fieldtripdefs

warning('this function is deprecated, please use the read_mri function instead');

% use the low-level reading function
[mri] = read_mri(filename);

