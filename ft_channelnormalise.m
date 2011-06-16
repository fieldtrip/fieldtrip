function [dataout] = ft_channelnormalise(cfg, data);

% FT_CHANNELNORMALISE normalises the input data to a mean of 0 and
% a standard deviation of 1.
%
% Use as
%   [dataout] = ft_channelnormalise(cfg, data)
%
% The configuration can contain
%   cfg.trials = 'all' or a selection given as a 1xN vector (default = 'all')
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
% Copyright (C) 2010, Jan-Mathijs Schoffelen

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

% set the defaults
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

% initialise some variables
nchan  = numel(data.label);
ntrl   = numel(data.trial);
datsum = zeros(nchan,1);
datssq = zeros(nchan,1);

% create output data, omitting sensor information
% FIXME this can be kept, provided the scaling is built in appropriately
dataout         = [];
dataout.label   = data.label;
if isfield(data, 'fsample'); dataout.fsample = data.fsample; end;
dataout.trial   = cell(1,ntrl);
dataout.time    = data.time;
if isfield(data, 'sampleinfo'),  dataout.sampleinfo  = data.sampleinfo;  end
if isfield(data, 'trialinfo'), dataout.trialinfo = data.trialinfo; end

% compute the mean and std
for k = 1:ntrl
    n(k,1) = size(data.trial{k},2);
    datsum = datsum + sum(data.trial{k},2);
    datssq = datssq + sum(data.trial{k}.^2,2);
end
datmean = datsum./sum(n);
datstd  = sqrt( (datssq - (datsum.^2)./sum(n))./sum(n)); %quick way to compute std from sum and sum-of-squared values

% demean and normalise
for k = 1:ntrl
  dataout.trial{k} = (data.trial{k}-datmean(:,ones(1,n(k))))./datstd(:,ones(1,n(k)));
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
dataout.cfg = cfg;

% convert back to input type if necessary
switch dtype 
    case 'timelock'
        dataout = ft_checkdata(dataout, 'datatype', 'timelock');
    otherwise
        % keep the output as it is
end

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', dataout); % use the variable name "data" in the output file
end
