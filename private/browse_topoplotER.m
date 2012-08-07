function browse_topoplotER(cfg, data)

% BROWSE_TOPOPLOTER is a simple helper function for FT_DATABROWSER and shows
% the average topography of the selected data.
%
% See also BROWSE_MOVIEPLOTER, BROWSE_TOPOPLOTER, BROWSE_MULTIPLOTER, BROWSE_TOPOPLOTVAR, BROWSE_SIMPLEFFT

% Copyright (C) 2009, Robert Oostenveld

% convert to an ERP
timelock = timelockanalysis([], data);

default             = [];
default.xlim        = [min(timelock.time) max(timelock.time)];
default.marker      = 'on';
default.interactive = 'no';
cfg = mergeconfig(cfg, default);

figure;
ft_topoplotER(cfg, timelock);
