function [interp] = ft_channelrepair(cfg, data);

% FT_CHANNELREPAIR repairs bad channels in MEG or EEG data by replacing them
% with the average of its neighbours. It cannot be used reliably to
% repair multiple bad channels that lie next to each other.
%
% Use as
%   [interp] = ft_channelrepair(cfg, data)
%
% The configuration can contain
%   cfg.badchannel     = cell-array, see FT_CHANNELSELECTION for details
%   cfg.neighbourdist  = default is 4 cm
%   cfg.trials         = 'all' or a selection given as a 1xN vector (default = 'all')
%
% Since a nearest neighbour average is used, the input should contain
% a gradiometer or electrode definition, i.e. data.grad or data.elec.
%
% See also FT_MEGREALIGN, FT_MEGPLANAR

% Undocumented local options:
%   cfg.inputfile        = one can specifiy preanalysed saved data as input
%   cfg.outputfile       = one can specify output as file to save to disk

% Copyright (C) 2004-2009, Robert Oostenveld

% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

% set the default configuration
if ~isfield(cfg, 'neighbourdist'), cfg.neighbourdist = 4;         end
if ~isfield(cfg, 'badchannel'),    cfg.badchannel = {};           end
if ~isfield(cfg, 'trials'),        cfg.trials = 'all';            end
if ~isfield(cfg, 'inputfile'),    cfg.inputfile = [];           end
if ~isfield(cfg, 'outputfile'),   cfg.outputfile = [];          end

hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    data = loadvar(cfg.inputfile, 'data');
  end
end

% check if the input data is valid for this function
data = checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data = ft_selectdata(data, 'rpt', cfg.trials);
end

% determine the type of data
iseeg = ft_senstype(data, 'eeg');
ismeg = ft_senstype(data, 'meg');

if iseeg
  sens = data.elec;
elseif ismeg
  sens = data.grad;
else
  error('the data should contain either an electrode or a gradiometer definition');
end

% get the selection of channels that are bad
cfg.badchannel = ft_channelselection(cfg.badchannel, data.label);
[goodchanlabels,goodchanindcs] = setdiff(data.label,cfg.badchannel);
[goodsenslabels,goodsensindcs] = intersect(sens.label,goodchanlabels);

Ntrials = length(data.trial);
Nchans = length(data.label);
Nsens  = length(sens.label);

repair = eye(Nchans,Nchans);
[badindx] = match_str(data.label, cfg.badchannel);

for k=badindx(:)'
  fprintf('repairing channel %s\n', data.label{k});
  repair(k,k) = 0;
  
  sensindx = match_str(sens.label, data.label{k});
  for l=goodsensindcs(:)'
    distance = norm(sens.pnt(l,:)-sens.pnt(sensindx,:));
    if distance<cfg.neighbourdist
      % include this channel as neighbour, weigh with inverse distance
      repair(k,l) = 1/distance;
      fprintf('  using neighbour %s\n', sens.label{l});
    end
  end
  
  % normalise the repair matrix to unit weight
  repair(k,:) = repair(k,:) ./ sum(repair(k,:));
end

% use sparse matrix to speed up computations
repair = sparse(repair);

% compute the repaired data for each trial
for i=1:Ntrials
  fprintf('repairing bad channels for trial %d\n', i);
  interp.trial{i} = repair * data.trial{i};
end

% store the realigned data in a new structure
interp.fsample = data.fsample;
interp.time    = data.time;
interp.label   = data.label;
if iseeg
  interp.elec  = sens;
else
  interp.grad  = sens;
end

if isfield(data, 'sampleinfo')
  interp.sampleinfo = data.sampleinfo;
end
if isfield(data, 'trialinfo')
  interp.trialinfo = data.trialinfo;
end

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% store the configuration of this function call, including that of the previous function call
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id$';
% remember the configuration details of the input data

% remember the configuration details of the input data
try cfg.previous = data.cfg;end

% remember the exact configuration details in the output
interp.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', interp); % use the variable name "data" in the output file
end
