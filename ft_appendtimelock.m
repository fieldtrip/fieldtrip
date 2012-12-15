function [timelock] = ft_appendtimelock(cfg, varargin)

% FT_APPENDTIMELOCK concatenates multiple timelock (ERP/ERF) data
% structures that have been processed seperately. If the input data
% structures contain different channels, it will be concatenated along the
% channel direction. If the channels are identical in the input data
% structures, the data will be concatenated along the repetition dimension.
%
% Use as
%   combined = ft_appendtimelock(cfg, timelock1, timelock2, ...)
%
% See also FT_TIMELOCKANALYSIS, FT_APPENDDATA, FT_APPENDFREQ, FT_APPENDSOURCE

% Copyright (C) 2011, Robert Oostenveld
%
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar    varargin
ft_preamble provenance varargin

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'timelock', 'feedback', 'yes', 'hassampleinfo', 'ifmakessense');
end

% set the defaults
cfg.channel = ft_getopt(cfg, 'channel', 'all');

% ensure that all inputs are sufficiently consistent
for i=1:length(varargin)
  if ~isequal(varargin{i}.time, varargin{1}.time)
    error('this function requires identical time axes for all input structures');
  end
end

% select the channels that are in every dataset
for i=1:length(varargin)
  cfg.channel = ft_channelselection(cfg.channel, varargin{i}.label);
end

% start with the initial output structure
timelock        = [];
timelock.time   = varargin{1}.time;
timelock.label  = cfg.channel;
timelock.dimord = 'rpt_chan_time';

nchan  = length(timelock.label);
ntime  = length(timelock.time);


if isfield(varargin{1}, 'trial')
  % these don't make sense when concatenating the avg
  hastrialinfo  = isfield(varargin{1}, 'triainfo');
  hassampleinfo = isfield(varargin{1}, 'sampleinfo');
  hascov        = isfield(varargin{1}, 'cov') && numel(size(varargin{1}.cov))==3;
  
  ntrial = zeros(size(varargin));
  for i=1:length(varargin)
    ntrial(i) = size(varargin{i}.trial, 1);
  end
  trialsel = cumsum([1 ntrial]);
  
  timelock.trial = zeros(sum(ntrial), nchan, ntime);
  if hastrialinfo,  timelock.trialinfo = zeros(sum(ntrial), size(varargin{1}.trialinfo,2)); end
  if hassampleinfo, timelock.sampleinfo = zeros(sum(ntrial), size(varargin{1}.sampleinfo,2)); end
  if hascov, timelock.cov = zeros(sum(ntrial), nchan, nchan); end
  
  for i=1:length(varargin)
    % copy the desired data into the output structure
    begtrial = trialsel(i);
    endtrial = trialsel(i+1)-1;
    chansel = match_str(varargin{i}.label, cfg.channel);
    timelock.trial(begtrial:endtrial,:,:) = varargin{i}.trial(:,chansel,:);
    if hastrialinfo,  timelock.trialinfo(begtrial:endtrial,:)   = varargin{i}.trialinfo(:,:); end
    if hassampleinfo, timelock.sampleinfo(begtrial:endtrial,:)  = varargin{i}.sampleinfo(:,:); end
    if hascov,        timelock.cov(begtrial:endtrial,:,:)       = varargin{i}.cov(:,chanselchansel); end
  end % for varargin
  
elseif isfield(varargin{1}, 'avg')
  hascov = isfield(varargin{1}, 'cov') && numel(size(varargin{1}.cov))==2;
  
  ntrial = numel(varargin);
  timelock.trial = zeros(ntrial, nchan, ntime);
  if hascov, timelock.cov = zeros(sum(ntrial),nchan,nchan); end
  
  for i=1:length(varargin)
    % copy the desired data into the output structure
    chansel = match_str(varargin{i}.label, cfg.channel);
    timelock.trial(i,:,:) = varargin{i}.avg(chansel,:);
    if hascov, timelock.cov(i,:,:) = varargin{i}.cov(chansel,chansel); end
  end % for varargin
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble previous   varargin
ft_postamble provenance timelock
ft_postamble history    timelock
ft_postamble savevar    timelock
ft_postamble debug
