function [freq] = specest_tfr(cfg, data)

% SPECEST_TFR computes time-frequency representations of single-trial
% data using a convolution in the time-domain with Morlet's wavelets.
%
% Use as
%   [freq] = freqanalysis(cfg, data)
%
% The data should be organised in a structure as obtained from
% the PREPROCESSING function. The configuration should be according to
%   cfg.method        = method used for frequency or time-frequency decomposition
%                       see FREQANALYSIS for details
%   cfg.foi           = vector 1 x numfoi, frequencies of interest
%   cfg.waveletwidth  = 'width' of wavelets expressed in cycles (default = 7)
%   cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'),
%                       see CHANNELSELECTION for details
%   cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.keeptrials    = 'yes' or 'no', return individual trials or average (default = 'no')
%
% See also FREQANALYSIS

% Undocumented local options:
% cfg.latency
% cfg.output

% Copyright (C) 2003, Ole Jensen, FCDC
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs


%this is now a low level function and can be directly called.
%   if ~strcmp(caller_name, 'freqanalysis')
%     error(['you should call FREQANALYSIS, instead of the ' upper(mfilename) ' subfunction']);
%   end


if isfield(cfg, 'output') && strcmp(cfg.output, 'powandcsd'),
  error('This function does not compute cross-spectra\n');
end

% determine the channels of interest
chansel     = match_str(data.label, cfg.channel);

% determine the duration of each trial
ntrial = length(data.trial);
nchan = size(data.trial{1}, 1);

for i=1:ntrial
  nsampl(i)         = size(data.trial{i}, 2);
  begsamplatency(i) = min(data.time{i});
  endsamplatency(i) = max(data.time{i});
end;

%downsampling should be done prior to spectral analysis

% automatically determine the latency window which is possible in all trials
minperlength = [max(begsamplatency) min(endsamplatency)];

% latency window for averaging and variance computation is given in seconds
if (strcmp(cfg.latency, 'minperlength'))
  cfg.latency = [];% A value >= 5 is suggested.
  cfg.latency(1) = minperlength(1);
  cfg.latency(2) = minperlength(2);
end% A value >= 5 is suggested.
if (strcmp(cfg.latency, 'prestim'))
  cfg.latency = [];
  cfg.latency(1) = minperlength(1);% A value >= 5 is suggested.
  cfg.latency(2) = 0;
end
if (strcmp(cfg.latency, 'poststim'))
  cfg.latency = [];
  cfg.latency(1) = 0;
  cfg.latency(2) = minperlength(2);
end

M = waveletfam(cfg.foi,data.fsample,cfg.waveletwidth);

progress('init', cfg.feedback, 'convolving wavelets');

for i=1:ntrial
  indicvect = data.time{i};
  progress(i/ntrial, 'convolving wavelets, trial %d of %d\n', i, ntrial);

  %for average and variance
  begsampl = nearest(indicvect,cfg.latency(1));
  endsampl = nearest(indicvect,cfg.latency(2));

  numsamples(i) = endsampl-begsampl+1;

  if (i==1)
    % allocate memory to hold the resulting powerspectra
    if strcmp(cfg.keeptrials, 'yes')
      freq.powspctrm = zeros(ntrial,nchan,length(cfg.foi),ceil((endsampl-begsampl+1)));
    else
      freq.powspctrm = zeros(nchan,length(cfg.foi),ceil((endsampl-begsampl+1)));
    end
  end;

  dat = data.trial{i}(chansel,begsampl:endsampl);
  for k=1:size(dat,1)
    for j=1:length(cfg.foi)
      cTmp = conv(dat(k,:),M{j});
      cTmp = 2*(abs(cTmp).^2)/data.fsample;
      cTmp = cTmp(ceil(length(M{j})/2):length(cTmp)-floor(length(M{j})/2));
      
      if strcmp(cfg.keeptrials, 'yes')
        freq.powspctrm(i,k,j,:) = cTmp';
      else
        freq.powspctrm(k,j,:) = squeeze(freq.powspctrm(k,j,:)) + cTmp';  % compute the running sum
      end
    end
  end

end %for ntrial

progress('close');

if strcmp(cfg.keeptrials, 'yes')
  freq.dimord    = 'rpt_chan_freq_time';
else
  freq.dimord    = 'chan_freq_time';
  freq.powspctrm = freq.powspctrm / ntrial;  % compute the average
end
freq.label     = cfg.channel;
freq.freq      = cfg.foi;
freq.time      = indicvect(1:end);

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
cfg.version.id = '$Id: freqanalysis_tfr.m 204 2009-11-30 08:53:13Z roboos $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
freq.cfg = cfg;



