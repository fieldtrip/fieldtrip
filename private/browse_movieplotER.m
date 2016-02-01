function browse_movieplotER(cfg, data)

% BROWSE_MOVIEPLOTER is a helper function for FT_DATABROWSER and makes a
% movie of the data that was selected. See ft_movieplotER for further details
% on the options that can be specified as cfg.selcfg in ft_databrowser.
%
% See also BROWSE_MOVIEPLOTER, BROWSE_TOPOPLOTER, BROWSE_MULTIPLOTER, BROWSE_TOPOPLOTVAR, BROWSE_SIMPLEFFT

% Copyright (c) 2009, Ingrid Nieuwenhuis

% convert to an ERP
timelock = ft_timelockanalysis([], data);

if isfield(cfg, 'framesfile')
  cfg.framesfile = [cfg.framesfile, '_s', num2str(data.cfg.trl(1)), 'tos', num2str(data.cfg.trl(2))];
end

figure('renderer','zbuffer');
ft_movieplotER(cfg, timelock);
