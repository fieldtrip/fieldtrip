function [freq] = freqbaseline(cfg, freq);

% FREQBASELINE performs baseline normalization for time-frequency data
%
% Use as
%    [freq] = freqbaseline(cfg, freq)
% where the freq data comes from FREQANALYSIS and the configuration
% should contain
%   cfg.baseline     = [begin end] (default = 'no')
%   cfg.baselinetype = 'absolute' 'relchange' 'relative' (default = 'absolute')
%
% See also FREQANALYSIS, TIMELOCKBASELINE

% Copyright (C) 2004-2006, Marcel Bastiaansen
% Copyright (C) 2005-2006, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

% check if the input data is valid for this function
freq = checkdata(freq, 'datatype', 'freq', 'feedback', 'yes');

% set the defaults
if ~isfield(cfg, 'baseline'),     cfg.baseline     = 'no';       end
if ~isfield(cfg, 'baselinetype'), cfg.baselinetype = 'absolute'; end % default is to use an absolute baseline

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
cfg.version.id = '$Id: freqbaseline.m,v 1.23 2009/01/20 13:01:31 sashae Exp $';
% remember the configuration details of the input data
try, cfg.previous = freq.cfg; end
% remember the exact configuration details in the output 
freq.cfg = cfg;

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

