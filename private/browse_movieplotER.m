function browse_movieplotER(cfg, data)

% this is a helper function for DATABROWSER which makes a movie of the data
% that was selected. See MOVIEPLOTER
%

% Copyright (C) 2009, Ingrid Nieuwenhuis

% Convert to an ERP
timelock = timelockanalysis([], data);

if isfield(cfg, 'framesfile')        
  cfg.framesfile = [cfg.framesfile, '_S', num2str(data.cfg.trl(1)), 'toS', num2str(data.cfg.trl(2))];
end

figure; 
movieplotER(cfg, timelock);