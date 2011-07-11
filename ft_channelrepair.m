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
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_MEGREALIGN, FT_MEGPLANAR

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

ft_defaults

% record start time and total processing time
ftFuncTimer = tic();
ftFuncClock = clock();

cfg = ft_checkconfig(cfg, 'trackconfig', 'on');

% set the default configuration
cfg = ft_checkconfig(cfg, 'required', {'neighbours'});

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

% store original datatype
dtype = ft_datatype(data);

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

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
connectivityMatrix = channelconnectivity(cfg, data);
connectivityMatrix = connectivityMatrix(:, goodchanindcs); % all chans x good chans

Ntrials = length(data.trial);
Nchans = length(data.label);
Nsens  = length(sens.label);

repair = eye(Nchans,Nchans);
[badindx] = match_str(data.label, cfg.badchannel);

for k=badindx'
    fprintf('repairing channel %s\n', data.label{k});
    repair(k,k) = 0;
    l = find(connectivityMatrix(k, :));
    % get bad channels out
    [a, b] = setdiff(data.label(l), data.label(badindx));
    l(~ismember(find(l), b)) = [];    
    % get corresponding ids for sens structure
    [a, b] = match_str(sens.label, data.label(l));
    goodsensindx = a(b);
    [a, b] = match_str(sens.label, data.label(k));
    badsensindx = a(b);
    fprintf('\tusing neighbour %s\n', sens.label{goodsensindx});
    distance = sqrt(sum((sens.pnt(goodsensindx,:) - repmat(sens.pnt(badsensindx, :), numel(goodsensindx), 1)).^2, 2));
    repair(k,l) = (1./distance);
    repair(k,l) = repair(k,l) ./ sum(repair(k,l));
end

% use sparse matrix to speed up computations
repair = sparse(repair);

% compute the repaired data for each trial
fprintf('\n');
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
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% store the configuration of this function call, including that of the previous function call
cfg.version.name = mfilename('fullpath');
cfg.version.id   = '$Id$';

% add information about the Matlab version used to the configuration
cfg.callinfo.matlab = version();
  
% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername();

% remember the configuration details of the input data
try cfg.previous = data.cfg; end

% remember the exact configuration details in the output
interp.cfg = cfg;

% convert back to input type if necessary
switch dtype 
    case 'timelock'
        interp = ft_checkdata(interp, 'datatype', 'timelock');
    otherwise
        % keep the output as it is
end

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', interp); % use the variable name "data" in the output file
end
