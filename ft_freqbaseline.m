function [freq] = ft_freqbaseline(cfg, freq)

% FT_FREQBASELINE performs baseline normalization for time-frequency data
%
% Use as
%    [freq] = ft_freqbaseline(cfg, freq)
% where the freq data comes from FT_FREQANALYSIS and the configuration
% should contain
%   cfg.baseline     = [begin end] (default = 'no')
%   cfg.baselinetype = 'absolute' 'relchange' 'relative' (default = 'absolute')
%
% See also FT_FREQANALYSIS, FT_TIMELOCKBASELINE, FT_FREQCOMPARISON

% Undocumented local options:
%   cfg.inputfile  = one can specifiy preanalysed saved data as input
%   cfg.outputfile = one can specify output as file to save to disk

% Copyright (C) 2004-2006, Marcel Bastiaansen
% Copyright (C) 2005-2006, Robert Oostenveld
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

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

% set the defaults
if ~isfield(cfg, 'baseline'),     cfg.baseline     = 'no';       end
if ~isfield(cfg, 'baselinetype'), cfg.baselinetype = 'absolute'; end % default is to use an absolute baseline
if ~isfield(cfg, 'inputfile'),    cfg.inputfile = [];            end
if ~isfield(cfg, 'outputfile'),   cfg.outputfile = [];           end

% load optional given inputfile as data
hasdata = (nargin>1);

if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    freq = loadvar(cfg.inputfile, 'data');
  end
end

% check if the input data is valid for this function
freq = ft_checkdata(freq, 'datatype', 'freq', 'feedback', 'yes');

% give a warning if the input is inconsistent
if ischar(cfg.baseline) && strcmp(cfg.baseline, 'no') && ~isempty(cfg.baselinetype)
  warning('no baseline correction done');
end

if ischar(cfg.baseline) && strcmp(cfg.baseline, 'yes')
  % default is to take the prestimulus interval
  cfg.baseline = [-inf 0];
elseif ischar(cfg.baseline) && strcmp(cfg.baseline, 'no')
  % nothing to do
  return
end

haspow = issubfield(freq, 'powspctrm');
hascoh = issubfield(freq, 'cohspctrm');

% we have to ensure that we don't end up with an inconsistent dataset
% remove cross-spectral densities since coherence cannot be computed any more
if isfield(freq, 'crsspctrm')
  freq = rmfield(freq, 'crsspctrm');
end
if isfield(freq, 'cohspctrmsem')
  freq = rmfield(freq, 'cohspctrmsem');
end
if isfield(freq, 'powspctrmsem')
  freq = rmfield(freq, 'powspctrmsem');
end

if strcmp(freq.dimord, 'chan_freq_time')
  % apply the desired method for the average, see the subfunctions below
  if strcmp(cfg.baselinetype, 'absolute')
    if haspow, freq.powspctrm = TFabschange(freq.time, freq.freq, freq.powspctrm, cfg.baseline); end
    if hascoh, freq.cohspctrm = TFabschange(freq.time, freq.freq, freq.cohspctrm, cfg.baseline); end
  elseif strcmp(cfg.baselinetype, 'relchange')
    if haspow, freq.powspctrm = TFrelchange(freq.time, freq.freq, freq.powspctrm, cfg.baseline); end
    if hascoh, freq.cohspctrm = TFrelchange(freq.time, freq.freq, freq.cohspctrm, cfg.baseline); end
  elseif strcmp(cfg.baselinetype, 'relative')
    if haspow, freq.powspctrm = TFrelative(freq.time, freq.freq, freq.powspctrm, cfg.baseline); end
    if hascoh, freq.cohspctrm = TFrelative(freq.time, freq.freq, freq.cohspctrm, cfg.baseline); end
    % elseif strcmp(cfg.baselinetype, 'zscore')
    %   freq.powspctrm = TFzscore(freq.time, freq.freq, freq.powspctrm,cfg.baseline);
  else
    error('unsupported method for baseline normalization');
  end
  
elseif strcmp(freq.dimord, 'rpt_chan_freq_time')
  % apply the desired method for each trial, see the subfunctions below
  if ~haspow || hascoh
    error('this only works for power, not for coherence');
  end
  
  Ntrial = size(freq.powspctrm,1);
  for i=1:Ntrial
    % Reshape freq.powspctrm into 3D matrix
    % This relies on dimord being 'rpt_chan_freq_time'
    tfdata = reshape(freq.powspctrm(i,:,:,:), ...
      size(freq.powspctrm,2), ...
      size(freq.powspctrm,3), ...
      size(freq.powspctrm,4));
    
    if strcmp(cfg.baselinetype, 'absolute'),
      freq.powspctrm(i,:,:,:) = TFabschange(freq.time, freq.freq, tfdata, cfg.baseline);
    elseif strcmp(cfg.baselinetype, 'relchange')
      freq.powspctrm(i,:,:,:) = TFrelchange(freq.time, freq.freq, tfdata, cfg.baseline);
    elseif strcmp(cfg.baselinetype, 'relative')
      freq.powspctrm(i,:,:,:) = TFrelative(freq.time, freq.freq, tfdata, cfg.baseline);
      % elseif strcmp(cfg.baselinetype, 'zscore')
      %   freq.powspctrm = TFzscore(freq.time, freq.freq, freq.powspctrm,cfg.baseline);
    else
      error('unsupported method for baseline normalization');
    end
  end
  
else
  error('unsupported data dimensions');
end

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id$';

% remember the configuration details of the input data
try, cfg.previous = freq.cfg; end

% remember the exact configuration details in the output
freq.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', freq); % use the variable name "data" in the output file
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [TFdata] = TFzscore(timeVec,freqVec,TFdata,TimeInt)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Compute relative change from baseline on a TFR representation as obtained from the framework software
% % NB fixed order of dimensions is assumed in the TFdata: channels, frequencies, timepoints
% tidx = find(timeVec >= TimeInt(1) & timeVec <= TimeInt(2));
% TFtmp = TFdata(:,:,tidx);
% for k=1:size(TFdata,2) % loop frequencies
%   for l=1:size(TFdata,1) % loop channels
%     TFbl   (l,k) = squeeze(mean(TFdata(l,k,tidx),3));     %compute average baseline power
%     TFblstd(l,k) = squeeze(std (TFdata(l,k,tidx),[], 3)); %compute standard deviation
%   end
% end
% for k=1:size(TFdata,2) % loop frequencies
%   for l=1:size(TFdata,1) % loop channels
%     TFdata(l,k,:) = ((TFdata(l,k,:) - TFbl(l,k)) / TFblstd(l,k));      % compute zscore
%   end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TFdata] = TFrelative(timeVec,freqVec,TFdata,TimeInt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute relative change from baseline on a TFR representation as obtained from the framework software
% NB fixed order of dimensions is assumed in the TFdata: channels, frequencies, timepoints
tidx = find(timeVec >= TimeInt(1) & timeVec <= TimeInt(2));

if length(size(TFdata))~=3,
  error('Time-frequency matrix should have three dimensions (chan,freq,time)');
end

for k=1:size(TFdata,2) % loop frequencies
  for l=1:size(TFdata,1) % loop channels
    TFbl(l,k) = nan_mean(TFdata(l,k,tidx),3);%compute average baseline power
    
    if TFbl(l,k) == 0,
      error('Average baseline power is zero');
    end
    
  end
end
for k=1:size(TFdata,2) % loop frequencies
  for l=1:size(TFdata,1) % loop channels
    TFdata(l,k,:) = TFdata(l,k,:) / TFbl(l,k);     % compute relative change (i.e. ratio)
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TFdata] = TFrelchange(timeVec,freqVec,TFdata,TimeInt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute relative change from baseline on a TFR representation as obtained from the framework software
% NB fixed order of dimensions is assumed in the TFdata: channels, frequencies, timepoints
tidx = find(timeVec >= TimeInt(1) & timeVec <= TimeInt(2));

if length(size(TFdata))~=3,
  error('Time-frequency matrix should have three dimensions (chan,freq,time)');
end

for k=1:size(TFdata,2) % loop frequencies
  for l=1:size(TFdata,1) % loop channels
    TFbl(l,k) = nan_mean(TFdata(l,k,tidx),3); %compute average baseline power
    
    if TFbl(l,k) == 0,
      error('Average baseline power is zero');
    end
    
  end
end

for k=1:size(TFdata,2) % loop frequencies
  for l=1:size(TFdata,1) % loop channels
    TFdata(l,k,:) = ((TFdata(l,k,:) - TFbl(l,k)) / TFbl(l,k)); % compute relative change
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TFdata] = TFabschange(timeVec,freqVec,TFdata,TimeInt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subtract baseline from a TFR representation as obtained from the framework software
% NB fixed order of dimensions is assumed in the TFdata: channels, frequencies, timepoints
tidx = find(timeVec >= TimeInt(1) & timeVec <= TimeInt(2));

if length(size(TFdata))~=3,
  error('Time-frequency matrix should have three dimensions (chan,freq,time)');
end

for k=1:size(TFdata,2) % loop frequencies
  for l=1:size(TFdata,1) % loop channels
    TFbl(l,k) = nan_mean(TFdata(l,k,tidx),3); %compute average baseline power
  end
end
for k=1:size(TFdata,2) % loop frequencies
  for l=1:size(TFdata,1) % loop channels
    TFdata(l,k,:) = TFdata(l,k,:) - TFbl(l,k);        % subtract baseline power
  end
end

