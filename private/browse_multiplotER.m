function browse_multiplotER(cfg, data)

% BROWSE_MULTIPLOTER is a simple helper function for FT_DATBROWSER and shows
% an interactive multiplot of the selected data.
%
% See also BROWSE_MOVIEPLOTER, BROWSE_TOPOPLOTER, BROWSE_MULTIPLOTER, BROWSE_TOPOPLOTVAR

% Copyright (C) 2009, Robert Oostenveld

% convert to an ERP
timelock = timelockanalysis([], data);

default             = [];
default.interactive = 'yes';
cfg = mergeconfig(cfg, default);

figure;
ft_multiplotER(cfg, timelock);
