function browse_movieplotER(cfg, data)

% browse_movieplotER is a helper function for ft_databrowser and makes a
% movie of the data that was selected. see ft_movieplotER for further details
% on the options that can be specified as cfg.selcfg in ft_databrowser.
%
% see also browse_movieplotER, browse_topoplotER, browse_multiplotER, browse_topoplotVAR

% Copyright (c) 2009, Ingrid Nieuwenhuis

% convert to an ERP
timelock = timelockanalysis([], data);

if isfield(cfg, 'framesfile')
  cfg.framesfile = [cfg.framesfile, '_s', num2str(data.cfg.trl(1)), 'tos', num2str(data.cfg.trl(2))];
end

figure;
ft_movieplotER(cfg, timelock);
