function [dat] = read_fcdc_data(filename, header, begsample, endsample, chanindx, continuous)

% this function is deprecated, please use the read_data function instead

% Copyright (C) 2003-2009, Robert Oostenveld, F.C. Donders Centre
%
% $Log: read_fcdc_data.m,v $
% Revision 1.45  2009/05/07 14:21:29  roboos
% deprecated the read_fcdc and write_fcdc functions, give warning and mention the correct function to be used
%

fieldtripdefs

warning('this function is deprecated, please use the read_data function instead');

% set the defaults for the optional input arguments
if nargin<5
  chanindx = [];
end
if nargin<6 || isempty(continuous)
  continuous = false;
end

% use the low-level reading function
[dat] = read_data(filename, 'header', header, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', ~continuous);

