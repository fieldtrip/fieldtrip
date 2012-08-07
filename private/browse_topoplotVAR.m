function browse_topoplotVAR(cfg, data)

% BROWSE_TOPOPLOTVAR is a simple helper function for FT_DATABROWSER that
% computes the variance of band-pass filtered data and makes a topographic
% plot. It serves to make a quick-and-dirty power topography.
%
% See also BROWSE_MOVIEPLOTER, BROWSE_TOPOPLOTER, BROWSE_MULTIPLOTER, BROWSE_TOPOPLOTVAR, BROWSE_SIMPLEFFT

% Copyright (C) 2009, Robert Oostenveld

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
ft_topoplotER(cfg, timelock);
