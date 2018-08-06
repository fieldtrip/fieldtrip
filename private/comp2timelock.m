function [timelock, cfg] = comp2timelock(cfg, comp)

% COMP2TIMELOCK transform the independent components into something
% on which the timelocked source reconstruction methods can
% perform their trick.

% Copyright (C) 2005, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% only convert, do not perform channel or component selection
timelock        = [];
timelock.avg    = comp.topo;
timelock.label  = comp.topolabel;
timelock.time   = 1:size(timelock.avg,2);
if isfield(comp, 'cfg'), timelock.cfg = comp.cfg; end
timelock.dimord = 'chan_time';

if isfield(comp, 'grad')
  timelock.grad = comp.grad;
end

if isfield(comp, 'elec')
  timelock.elec = comp.elec;
end

