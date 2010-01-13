function [interp] = channelrepair(cfg, data);

% CHANNELREPAIR repairs bad channels in MEG or EEG data by replacing them
% with the average of its neighbours. It cannot be used reliably to
% repair multiple bad channels that ly next to each other.
%
% Use as
%   [interp] = channelrepair(cfg, data)
%
% The configuration can contain
%   cfg.badchannel     = cell-array, see CHANNELSELECTION for details
%   cfg.neighbourdist  = default is 4 cm 
%   cfg.trials         = 'all' or a selection given as a 1xN vector (default = 'all')
%
% Since a nearest neighbour average is used, the input should contain
% a gradiometer or electrode definition, i.e. data.grad or data.elec.
%
% See also MEGINTERPOLATE

% Copyright (C) 2004-2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

% check if the input data is valid for this function
data = checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

% set the default configuration 
if ~isfield(cfg, 'neighbourdist'), cfg.neighbourdist = 4;         end
if ~isfield(cfg, 'badchannel'),    cfg.badchannel = {};           end
if ~isfield(cfg, 'trials'),        cfg.trials = 'all';            end

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data = selectdata(data, 'rpt', cfg.trials);  
  % update the trial definition (trl)
  if isfield(data, 'cfg') % try to locate the trl in the nested configuration
    cfg.trl    = findcfg(data.cfg, 'trl');
    cfg.trlold = findcfg(data.cfg, 'trlold');
  end
end

% determine the type of data
iseeg = senstype(data, 'eeg');
ismeg = senstype(data, 'meg');

if iseeg
  sens = data.elec;
elseif ismeg
  sens = data.grad;
else
  error('the data should contain either an electrode or a gradiometer definition');
end

% get the selection of channels that are bad
cfg.badchannel = channelselection(cfg.badchannel, data.label);
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
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output 
interp.cfg = cfg;

