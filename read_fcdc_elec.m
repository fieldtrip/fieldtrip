function [sens] = read_fcdc_elec(filename)

% this function is deprecated, please use the read_sens function instead

% Copyright (C) 2005-200, Robert Oostenveld
%
% $Log: read_fcdc_elec.m,v $
% Revision 1.12  2009/05/07 14:21:29  roboos
% deprecated the read_fcdc and write_fcdc functions, give warning and mention the correct function to be used
%

fieldtripdefs

warning('this function is deprecated, please use the read_sens function instead');

% use the low-level reading function
[sens] = read_sens(filename);

