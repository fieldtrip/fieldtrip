function [dat] = read_fcdc_data(filename, header, begsample, endsample, chanindx, continuous)

% this function is deprecated, please use the read_data function instead

% Copyright (C) 2003-2009, Robert Oostenveld, F.C. Donders Centre
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

ft_defaults

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

