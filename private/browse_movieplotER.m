function browse_movieploter(cfg, data)

% browse_movieploter is a helper function for ft_databrowser and makes a
% movie of the data that was selected. see ft_movieploter for further details
% on the options that can be specified as cfg.selcfg in ft_databrowser.
%
% see also browse_movieploter, browse_topoploter, browse_multiploter, browse_topoplotvar

% Copyright (c) 2009, ingrid nieuwenhuis

% convert to an ERP
timelock = timelockanalysis([], data);

if isfield(cfg, 'framesfile')
  cfg.framesfile = [cfg.framesfile, '_s', num2str(data.cfg.trl(1)), 'tos', num2str(data.cfg.trl(2))];
end

figure;
ft_movieploter(cfg, timelock);