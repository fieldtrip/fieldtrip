function [timelock, cfg] = comp2timelock(cfg, comp);

% COMP2TIMELOCK transform the independent components into something
% on which the timelocked source reconstruction methods can
% perform their trick.

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: comp2timelock.m,v $
% Revision 1.2  2006/05/10 08:19:45  roboos
% added dimord to the output
%
% Revision 1.1  2005/10/14 15:50:08  roboos
% new implementation, used by dipolefitting in case of frequency or ICA data
%

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

