function [data] = megrepair(cfg, data)

% MEGREPAIR is deprecated, please use CHANNELREPAIR

% Copyright (C) 2003-2009, Robert Oostenveld, F.C. Donders Centre
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

warning('MEGREPAIR is deprecated, please use CHANNELREPAIR');
data = channelrepair(cfg, data);

