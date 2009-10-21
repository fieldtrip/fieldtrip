function [data] = megrepair(cfg, data)

% MEGREPAIR is deprecated, please use CHANNELREPAIR

% Copyright (C) 2003-2009, Robert Oostenveld, F.C. Donders Centre
%
% $Log: megrepair.m,v $
% Revision 1.21  2009/02/04 09:33:49  roboos
% renamed megrepair to channelrepair, to make clear that it is not a MEG specific function
%

warning('MEGREPAIR is deprecated, please use CHANNELREPAIR');
data = channelrepair(cfg, data);

