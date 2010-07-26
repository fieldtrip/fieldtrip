function [freq] =ft_freqanalysis(cfg, data, flag);

% FT_FREQANALYSIS performs frequency and time-frequency analysis
% on time series data over multiple trials
%
% Use as
%   [freq] = ft_freqanalysis(cfg, data)
%
% The input data should be organised in a structure as obtained from
% the FT_PREPROCESSING function. The configuration depends on the type
% of computation that you want to perform.
%
% The configuration should contain:
%   cfg.method     = different methods of calculating the spectra
%                    'mtmfft', analyses an entire spectrum for the entire data
%                     length, implements multitaper frequency transformation
%                    'mtmconvol', implements multitaper time-frequency transformation
%                     based on multiplication in the frequency domain
%                    'mtmwelch', performs frequency analysis using Welch's averaged
%                     modified periodogram method of spectral estimation
%                    'wltconvol', implements wavelet time frequency transformation
%                     (using Morlet wavelets) based on multiplication in the frequency domain
%                    'tfr', implements wavelet time frequency transformation
%                     (using Morlet wavelets) based on convolution in the time domain
%
%
% The other cfg options depend on the method that you select. You should
% read the help of the respective subfunction FT_FREQANALYSIS_XXX for the
% corresponding parameter options and for a detailed explanation of each method.
%
% See also FT_FREQANALYSIS_MTMFFT, FT_FREQANALYSIS_MTMCONVOL, FT_FREQANALYSIS_MTMWELCH
% FT_FREQANALYSIS_WLTCONVOL, FT_FREQANALYSIS_TFR

% Undocumented local options:
% cfg.channel
% cfg.channelcmb
% cfg.inputfile  = one can specifiy preanalysed saved data as input
% cfg.outputfile = one can specify output as file to save to disk

% Copyright (C) 2003-2006, F.C. Donders Centre, Pascal Fries
% Copyright (C) 2004-2006, F.C. Donders Centre, Markus Siegel
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

%allow for both the new and old implementation to be changed with a flag
%input

% defaults for optional input/ouputfile
if ~isfield(cfg, 'inputfile'),  cfg.inputfile                   = [];    end
if ~isfield(cfg, 'outputfile'), cfg.outputfile                  = [];    end

% load optional given inputfile as data
hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    data = loadvar(cfg.inputfile, 'data');
  end
end

if nargin < 3
  flag = 0;
end

% check if the input data is valid for this function
data = checkdata(data, 'datatype', {'raw', 'comp', 'mvar'}, 'feedback', 'yes', 'hasoffset', 'yes', 'hastrialdef', 'yes');

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'renamed',     {'label', 'channel'});
cfg = checkconfig(cfg, 'renamed',     {'sgn',   'channel'});
cfg = checkconfig(cfg, 'renamed',     {'labelcmb', 'channelcmb'});
cfg = checkconfig(cfg, 'renamed',     {'sgncmb',   'channelcmb'});
cfg = checkconfig(cfg, 'required',    {'method'});
cfg = checkconfig(cfg, 'renamedval',  {'method', 'fft',    'mtmfft'});
cfg = checkconfig(cfg, 'renamedval',  {'method', 'convol', 'mtmconvol'});

% select trials of interest
if ~isfield(cfg, 'trials'),   cfg.trials = 'all';  end % set the default
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data = selectdata(data, 'rpt', cfg.trials);
end

if ~flag
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % HERE THE OLD IMPLEMENTATION STARTS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [freq] = feval(sprintf('ft_freqanalysis_%s',lower(cfg.method)), cfg, data);
  
else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % HERE THE NEW IMPLEMENTATION STARTS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % set all the defaults
  if ~isfield(cfg, 'padding'), cfg.padding = [];   end
  if ~isfield(cfg, 'output'),  cfg.output = 'pow'; end
  if ~isfield(cfg, 'taper'),   cfg.taper =  'dpss';     end
  if ~isfield(cfg, 'method'), error('you must specify a method in cfg.method'); end
  if isequal(cfg.taper, 'dpss') && not(isfield(cfg, 'tapsmofrq'))
    error('you must specify a smoothing parameter with taper = dpss');
  end
  if ~isfield(cfg, 'keeptapers'),    cfg.keeptapers = 'no';         end
  if ~isfield(cfg, 'keeptrials'),    cfg.keeptrials = 'no';         end
  if ~isfield(cfg, 'calcdof'),       cfg.calcdof    = 'no';         end
  
  if ~isfield(cfg, 'pad'),           cfg.pad        = 'maxperlen';  end
  if ~isfield(cfg, 'channel'),       cfg.channel    = 'all';        end
  if ~isfield(cfg, 'precision'),     cfg.precision  = 'double';     end
  if ~isfield(cfg, 'output'),        cfg.output     = 'powandcsd';  end
  if strcmp(cfg.output, 'fourier'),
    cfg.keeptrials = 'yes';
    cfg.keeptapers = 'yes';
  end
  if ~isfield(cfg, 'foi'),        cfg.foi    = [];    end
  if ~isfield(cfg, 'foilim'),     cfg.foilim = [];    end
  
 
  % set flags for keeping trials and/or tapers
  if strcmp(cfg.keeptrials,'no') &&  strcmp(cfg.keeptapers,'no')
    keeprpt = 1;
  elseif strcmp(cfg.keeptrials,'yes') &&  strcmp(cfg.keeptapers,'no')
    keeprpt = 2;
  elseif strcmp(cfg.keeptrials,'no') &&  strcmp(cfg.keeptapers,'yes')
    error('There is currently no support for keeping tapers WITHOUT KEEPING TRIALS.');
  elseif strcmp(cfg.keeptrials,'yes') &&  strcmp(cfg.keeptapers,'yes')
    keeprpt = 4;
  end
  
  % Set flags for output
  if strcmp(cfg.output,'pow')
    powflg = 1;
    csdflg = 0;
    fftflg = 0;
  elseif strcmp(cfg.output,'powandcsd')
    powflg = 1;
    csdflg = 1;
    fftflg = 0;
  elseif strcmp(cfg.output,'fourier')
    powflg = 0;
    csdflg = 0;
    fftflg = 1;
  else
    error('Unrecognized output required');
  end
  
  % prepare channel(cmb)
  if ~isfield(cfg, 'channelcmb') && csdflg
    %set the default for the channelcombination
    cfg.channelcmb = {'all' 'all'};
  elseif isfield(cfg, 'channelcmb') && ~csdflg
    % no cross-spectrum needs to be computed, hence remove the combinations from cfg
    cfg = rmfield(cfg, 'channelcmb');
  end
  
  % ensure that channelselection and selection of channelcombinations is
  % perfomed consistently
  cfg.channel = ft_channelselection(cfg.channel, data.label);
  if isfield(cfg, 'channelcmb')
    cfg.channelcmb = ft_channelcombination(cfg.channelcmb, data.label);
  end
  
  % determine the corresponding indices of all channels
  chanind    = match_str(data.label, cfg.channel);
  nchan      = size(chanind,1);
  if csdflg
    % determine the corresponding indices of all channel combinations
    chancmbind = zeros(size(cfg.channelcmb));
    for k=1:size(cfg.channelcmb,1)
      chancmbind(k,1) = strmatch(cfg.channelcmb(k,1), data.label, 'exact');
      chancmbind(k,2) = strmatch(cfg.channelcmb(k,2), data.label, 'exact');
    end
    nchancmb   = size(chancmbind,1);
    chanind    = unique([chanind(:); chancmbind(:)]);
    nchan      = length(chanind);
    cutdatindcmb = zeros(size(chancmbind));
    for ichan = 1:nchan
      cutdatindcmb(find(chancmbind == chanind(ichan))) = ichan;
    end
  end
  
  % determine trail characteristics
  ntrials = size(data.trial,2);
  trllength = zeros(1, ntrials);
  for itrial = 1:ntrials
    trllength(itrial) = size(data.trial{itrial}, 2);
  end
  if strcmp(cfg.padding, 'maxperlen')
    padding = max(trllength);
    cfg.padding = padding/data.fsample;
  else
    padding = cfg.padding*data.fsample;
    cfg.padding = padding;
    if padding<max(trllength)
      error('the specified padding is too short');
    end
  end
  
  % foilim 'backwards compatibility'
  if ~isempty(cfg.foilim) && isempty(cfg.foilim)
    error('use either cfg.foi or cfg.foilim')
  elseif ~isempty(cfg.foilim)
    cfg.foi = cfg.foilim(1):data.fsample/(cfg.padding):cfg.foilim(2); % get the full foi in the current foilim and set it too be used as foilim
  end
  
  % options that don't change over trials
  options = {'pad', cfg.padding, 'taper', cfg.taper, 'freqoi', cfg.foi, 'tapsmofrq', cfg.tapsmofrq};
  
    
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Main loop over trials, inside fourierspectra are obtained and transformed into the appropriate outputs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % this is done on trial basis to save memory
  for itrial = 1:ntrials
    disp(['processing trial ' num2str(itrial) ': ' num2str(size(data.trial{itrial},2)) ' samples']);
    
    dat = data.trial{itrial}; % chansel has already been performed
    time = data.time{itrial};
    
    % Perform memory allocation and several other things that require at least 1 specest call
    if itrial == 1
      % minimal specest call, bookkeeping and such is in a second switch below here
      foi = [];
      toi = []; % declare
      switch cfg.method
        case 'mtmconvol'
          [spectrum,foi,toi] = specest_mtmconvol(dat, time, 'timeoi', cfg.toi, options{:},'timwin',cfg.t_ftimwin);
        case 'mtmfft'
          [spectrum,foi] = specest_mtmfft(dat, time, options{:});
        otherwise
          error('method %s is unknown', cfg.method);
      end % switch
      if ~isempty(foi)
        nfoi = numel(foi);
      end
      if ~isempty(toi)
        ntoi = numel(toi);
      end
      ntap = size(spectrum,1);
      
      % allocate memory to output variables
      if keeprpt == 1 % cfg.keeptrials,'no' &&  cfg.keeptapers,'no'
        if powflg, powspctrm     = zeros(nchan,nfoi,ntoi,cfg.precision);             end
        if csdflg, crsspctrm     = complex(zeros(nchancmb,nfoi,ntoi,cfg.precision)); end
        if fftflg, fourierspctrm = complex(zeros(nchan,nfoi,ntoi,cfg.precision));    end
        dimord    = 'chan_freq_time';
      elseif keeprpt == 2 % cfg.keeptrials,'yes' &&  cfg.keeptapers,'no'
        if powflg, powspctrm     = zeros(ntrials,nchan,nfoi,ntoi,cfg.precision);             end
        if csdflg, crsspctrm     = complex(zeros(ntrials,nchancmb,nfoi,ntoi,cfg.precision)); end
        if fftflg, fourierspctrm = complex(zeros(ntrials,nchan,nfoi,ntoi,cfg.precision));    end
        dimord    = 'rpt_chan_freq_time';
      elseif keeprpt == 4 % cfg.keeptrials,'yes' &&  cfg.keeptapers,'yes'
        % FIXME this works only if all frequencies have the same number of tapers
        if powflg, powspctrm     = zeros(ntrials*ntap,nchan,nfoi,ntoi,cfg.precision);             end
        if csdflg, crsspctrm     = complex(zeros(ntrials*ntap,nchancmb,nfoi,ntoi,cfg.precision)); end
        if fftflg, fourierspctrm = complex(zeros(ntrials*ntap,nchan,nfoi,ntoi,cfg.precision));    end
        dimord    = 'rpttap_chan_freq_time';
      end
    end % if itrial == 1
    
    
    
    
    
    % do the spectral decompisition of this trial and put it into output variables
    switch cfg.method
      
      case 'mtmconvol'
        if itrial ~= 1
          [spectrum,foi,toi] = specest_mtmconvol(dat, time, 'timeoi', cfg.toi, options{:},'timwin',cfg.t_ftimwin);
          nfoi = numel(foi);
          ntoi = numel(toi);
          ntap = size(spectrum,1);
        end
        
        % get output in correct format
        % for now, there is a lot of redundancy, as each method has it's own case statement
        % when fully implemented, this can be cut down, perhaps in a separate switch, or perhaps as a time and a non-time if-loop
        foinumsmp = cfg.t_ftimwin .* data.fsample;
        foinumsmp = reshape(repmat(foinumsmp,[1, ntap, nchan, ntoi]),[ntap, nchan, nfoi, ntoi]);
        if powflg
          powdum = 2.* abs(spectrum) .^ 2 ./ foinumsmp;
          %           if strcmp(cfg.taper, 'sine') % NOT YET RE-IMPLEMENTED
          %             powdum = powdum .* (1 - (((taplop - 1) ./ ntap(foilop)) .^ 2));
          %           end
          if keeprpt == 1 % cfg.keeptrials,'no' &&  cfg.keeptapers,'no'
            powspctrm = powspctrm + (reshape(mean(powdum,1),[nchan nfoi ntoi]) ./ ntrials);
          elseif keeprpt == 2  % cfg.keeptrials,'yes' &&  cfg.keeptapers,'no'
            powspctrm(itrial,:,:,:) = reshape(mean(powdum,1),[nchan nfoi ntoi]);
          elseif keeprpt == 4 % cfg.keeptrials,'yes' &&  cfg.keeptapers,'yes'
            rptind = reshape(1:ntrials * ntap,[ntap ntrials]);
            powspctrm(rptind(:,itrial),:,:,:) = reshape(powdum,[ntap nchan nfoi ntoi]);
          end
        end
        if fftflg
          fourierdum = spectrum .* sqrt(2 ./ foinumsmp); %cf Numercial Receipes 13.4.9
          if keeprpt == 1
            fourierspctrm = fourierspctrm + (reshape(mean(fourierdum,1),[nchan nfoi ntoi]) ./ ntrials);
          elseif keeprpt == 2
            fourierspctrm(itrial,:,:,:) = reshape(mean(fourierdum,1), [nchan nfoi ntoi]);
          elseif keeprpt == 4
            rptind = reshape(1:ntrials * ntap,[ntap ntrials]);
            fourierspctrm(rptind(:,itrial),:,:,:) = reshape(fourierdum,[ntap nchan nfoi ntoi]);
          end
        end
        if csdflg
          csddum = 2.* (spectrum(:,cutdatindcmb(:,1),:,:) .* conj(spectrum(:,cutdatindcmb(:,2),:,:))) ./ foinumsmp;
          if keeprpt == 1
            crsspctrm = crsspctrm + (reshape(mean(csddum,1),[nchancmb nfoi ntoi]) ./ ntrials);
          elseif keeprpt == 2
            crsspctrm(itrial,:,:,:) = reshape(mean(csddum,1), [nchancmb nfoi ntoi]);
          elseif keeprpt == 4
            rptind = reshape(1:ntrials * ntap,[ntap ntrials]);
            crsspctrm(rptind(:,itrial),:,:,:) = reshape(csddum,[ntap nchancmb nfoi ntoi]);
          end
        end
        
        
      case 'mtmfft'
        if itrial ~= 1
          [spectrum,ntap,foi] = specest_mtmfft(dat, time, options{:});
        end
        
        %       case 'wltconvol' % not testest yet
        %         [spectrum,foi,toi] = specest_wltconvol(dat, time, options{:});
        
        
      otherwise
        error('method %s is unknown', cfg.method);
    end % switch
    
    
    
    
    
    
    %     % now get the output in the correct format
    %     if strcmp(cfg.output, 'pow')
    %       freq.dimord = 'rpt_chan_freq';
    %       freq.freq = foi;
    %       freq.powspctrm = 2 .* (fourierspctrm .* conj(fourierspctrm)) ./ numsmp; % cf Numercial Receipes 13.4.9
    %       freq.powspctrm = reshape(freq.powspctrm, size(freq.powspctrm,1)*size(freq.powspctrm,2),size(freq.powspctrm,3),size(freq.powspctrm,4));
    %     elseif strcmp(cfg.output, 'fourier')
    %       freq.dimord = 'rpt_chan_freq';
    %       freq.freq = foi;
    %       freq.fourierspctrm = (fourierspctrm) .* sqrt(2 ./ numsmp); % cf Numercial Receipes 13.4.9
    %       freq.fourierspctrm = reshape(freq.fourierspctrm, size(freq.fourierspctrm,1)*size(freq.fourierspctrm,2),size(freq.fourierspctrm,3),size(freq.fourierspctrm,4));
    %     elseif strcmp(cfg.output, 'csd')
    %       freq.dimord = 'rpt_chan_freq';
    %       freq.freq = foi;
    %       freq.cumtapcnt = size(spectrum,1)*zeros(ntrials,1);
    %       freq.fourierspctrm = (fourierspctrm) .* sqrt(2 ./ numsmp); % cf Numercial Receipes 13.4.9
    %       freq.fourierspctrm = reshape(freq.fourierspctrm, size(freq.fourierspctrm,1)*size(freq.fourierspctrm,2),size(freq.fourierspctrm,3),size(freq.fourierspctrm,4));
    %       freq = fixcsd(freq, 'full', []);
    %     else
    %       error('output is not recognized',cfg.output);
    %     end
    
  end % for ntrials
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% END: Main loop over trials
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  
  % set output variables
  freq = [];
  freq.label = data.label;
  freq.dimord = dimord;
  freq.freq   = foi;
  if exist('toi','var')
    freq.time = toi;
  end
  if powflg
    freq.powspctrm = powspctrm;
  end
  if fftflg
    freq.fourierspctrm = fourierspctrm;
  end
  if csdflg
    freq.crsspctrm = crsspctrm;
  end
  
  
  
  
  
  % accessing this field here is needed for the configuration tracking
  % by accessing it once, it will not be removed from the output cfg
  cfg.outputfile;
  
  % get the output cfg
  cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');
  
  % add information about the version of this function to the configuration
  try
    % get the full name of the function
    cfg.version.name = mfilename('fullpath');
  catch
    % required for compatibility with Matlab versions prior to release 13 (6.5)
    [st, i1] = dbstack;
    cfg.version.name = st(i1);
  end
  cfg.version.id = '$Id$';
  % remember the configuration details of the input data
  try, cfg.previous = data.cfg; end
  % remember the exact configuration details in the output
  freq.cfg = cfg;
  
end % IF OLD OR NEW IMPLEMENTATION









% copy the trial specific information into the output
if isfield(data, 'trialinfo'),
  freq.trialinfo = data.trialinfo;
end

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', freq); % use the variable name "data" in the output file
end






