function realtime_packettimer(cfg)

% REALTIME_PACKETTIMER can be used to time the rate at which data can be processed
%
% Use as
%   realtime_packettimer(cfg)
% with the following configuration options
%   cfg.bcifun    = processing of the data (default = @bcifun_timer)
%   cfg.npackets  = the number of packets shown in one plot (default=1000)
%                     after reaching the end
%   cfg.saveplot  = if path is specified, first plot is saved (default=[]);
%   cfg.rellim = y limits of subplot 1 (default = [-100 100])
%
% SEE ALSO:
%   FT_REALTIME_PROCESS

% TO DO:
%   jitter in het binnenhalen van de data; scatterplot!
%   triggers sturen en herhalen (loop closen)
%   tijd schatten waarin matlab nog kan processen

% Copyright (C) 2009, Marcel van Gerven
%
% $Id$

if ~isfield(cfg,'bcifun'),    cfg.bcifun = @bcifun_timer; end
if ~isfield(cfg,'npackets'),  cfg.npackets = 10^2; end
if ~isfield(cfg,'rellim'),    cfg.rellim = [-1 1]; end
if ~isfield(cfg,'saveplot'),  cfg.saveplot= []; end

close all;
f1 = figure();
set(f1,'units','normalized','outerposition',[0 0 1 1]);

% reset persistent variables
cfg.bcifun();

try
  realtime_process(cfg);
catch
  fprintf('%s\n',lasterr);
  close;
end

