function browse_multiplotER(cfg, data)

% BROWSE_MULTIPLOTER is a simple helper function for FT_DATABROWSER and shows
% an interactive multiplot of the selected data.
%
% See also BROWSE_MOVIEPLOTER, BROWSE_TOPOPLOTER, BROWSE_MULTIPLOTER, BROWSE_TOPOPLOTVAR, BROWSE_SIMPLEFFT

% Copyright (C) 2009, Robert Oostenveld

% convert to an ERP
timelock = ft_timelockanalysis([], data);

default             = [];
default.interactive = 'yes';
default.axes        = 'no';
default.baseline    = 'yes';
cfg = mergestruct(cfg, default);

figure;
ft_multiplotER(cfg, timelock);
