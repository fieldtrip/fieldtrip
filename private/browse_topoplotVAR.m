function browse_topoplotVAR(cfg, data)

% this is a simple helper function for DATABROWSER

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: browse_topoplotVAR.m,v $
% Revision 1.1  2009/10/19 14:19:21  roboos
% first version, to work with databrowser
%

% compute the variance, i.e. the broad-band power
timelock        = [];
timelock.label  = data.label;
timelock.time   = mean(data.time{1});
timelock.avg    = sum(preproc_baselinecorrect(data.trial{1}).^2, 2);
timelock.dimord = 'chan_time';

default             = [];
default.markers     = 'labels';
default.interactive = 'no';
cfg = mergeconfig(cfg, default);

figure;
topoplotER(cfg, timelock);
