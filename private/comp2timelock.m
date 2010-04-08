function [timelock, cfg] = comp2timelock(cfg, comp);

% COMP2TIMELOCK transform the independent components into something
% on which the timelocked source reconstruction methods can
% perform their trick.

% Copyright (C) 2005, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% only convert, do not perform channel or component selection
timelock        = [];
timelock.avg    = comp.topo;
timelock.label  = comp.topolabel;
timelock.time   = 1:size(timelock.avg,2);
timelock.cfg    = comp.cfg;
timelock.dimord = 'chan_time';

if isfield(comp, 'grad')
  timelock.grad = comp.grad;
end

if isfield(comp, 'elec')
  timelock.elec = comp.elec;
end

