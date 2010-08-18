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
% cfg.correctt_ftimwin (set to yes to try to determine new t_ftimwins based on correct cfg.foi)
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
  if ~isfield(cfg, 'pad'),              cfg.pad              = 'maxperlen';  end
  if ~isfield(cfg, 'output'),           cfg.output           = 'pow';        end
  if ~isfield(cfg, 'taper'),            cfg.taper            =  'dpss';      end
  if ~isfield(cfg, 'method'), error('you must specify a method in cfg.method'); end
  if isequal(cfg.taper, 'dpss') && not(isfield(cfg, 'tapsmofrq'))
    error('you must specify a smoothing parameter with taper = dpss');
  end
  if ~isfield(cfg, 'keeptapers'),       cfg.keeptapers       = 'no';         end
  if ~isfield(cfg, 'keeptrials'),       cfg.keeptrials       = 'no';         end
  if ~isfield(cfg, 'calcdof'),          cfg.calcdof          = 'no';         end
  
  if ~isfield(cfg, 'pad'),              cfg.pad              = 'maxperlen';  end
  if ~isfield(cfg, 'channel'),          cfg.channel          = 'all';        end
  if ~isfield(cfg, 'precision'),        cfg.precision        = 'double';     end
  if ~isfield(cfg, 'output'),           cfg.output           = 'powandcsd';  end
  if strcmp(cfg.output, 'fourier'),
    cfg.keeptrials = 'yes';
    cfg.keeptapers = 'yes';
  end
  if ~isfield(cfg, 'foi'),              cfg.foi              = [];           end
  if ~isfield(cfg, 'foilim'),           cfg.foilim           = [];           end
  if ~isfield(cfg, 'correctt_ftimwin'), cfg.correctt_ftimwin = 'no';         end
  
  
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
  ntrials = numel(data.trial);
  trllength = zeros(1, ntrials);
  for itrial = 1:ntrials
    trllength(itrial) = size(data.trial{itrial}, 2);
  end
  if strcmp(cfg.pad, 'maxperlen')
    padding = max(trllength);
    cfg.pad = padding/data.fsample;
  else
    padding = cfg.pad*data.fsample;
    cfg.pad = padding;
    if padding<max(trllength)
      error('the specified padding is too short');
    end
  end
  
  % correct foi and implement foilim 'backwards compatibility'
  if ~isempty(cfg.foi) && ~isempty(cfg.foilim)
    error('use either cfg.foi or cfg.foilim')
  elseif ~isempty(cfg.foilim)
    % get the full foi in the current foilim and set it too be used as foilim
    fboilim = round(cfg.foilim ./ (data.fsample ./ (cfg.pad*data.fsample))) + 1;
    fboi    = fboilim(1):1:fboilim(2);
    cfg.foi = (fboi-1) ./ cfg.pad;
  else
    % correct foi if foilim was empty and try to correct t_ftimwin (by detecting whether there is a constant factor between foi and t_ftimwin: cyclenum)
    oldfoi = cfg.foi;
    fboi   = round(cfg.foi ./ (data.fsample ./ (cfg.pad*data.fsample))) + 1;
    cfg.foi    = (fboi-1) ./ cfg.pad; % boi - 1 because 0 Hz is included in fourier output
    if strcmp(cfg.correctt_ftimwin,'yes')
      cyclenum = oldfoi .* cfg.t_ftimwin;
      cfg.t_ftimwin = cyclenum ./ cfg.foi;
    end
  end
  
  
  % correct t_ftimwin to proper foi
  
  % tapsmofrq compatibility between functions (make it into a vector if it's not)
  if isfield(cfg,'tapsmofrq')
    if strcmp(cfg.method,'mtmconvol') && length(cfg.tapsmofrq) == 1 && length(cfg.foi) ~= 1
      cfg.tapsmofrq = ones(length(cfg.foi),1) * cfg.tapsmofrq;
    elseif strcmp(cfg.method,'mtmfft') && length(cfg.tapsmofrq) ~= 1
      cfg.tapsmofrq = cfg.tapsmofrq(1);
    end
  end
  
  
  % options that don't change over trials
  if isfield(cfg,'tapsmofrq')
    options = {'pad', cfg.pad, 'taper', cfg.taper, 'freqoi', cfg.foi, 'tapsmofrq', cfg.tapsmofrq};
  else
    options = {'pad', cfg.pad, 'taper', cfg.taper, 'freqoi', cfg.foi};
  end
  
  
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
      switch cfg.method
        case 'mtmconvol'
          [spectrum,ntaper,foi,toi] = specest_mtmconvol(dat, time, 'timeoi', cfg.toi, 'timwin', cfg.t_ftimwin, options{:});
          nfoi = numel(foi);
          ntoi = numel(toi);
          ntap = size(spectrum,1);
          if strcmp(cfg.calcdof,'yes')
            dof = zeros(ntrials,nfoi,ntoi);
          end
          
          
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
            % FIXME this works only if all frequencies have the same number of tapers (ancient fixme)
            if powflg, powspctrm     = zeros(ntrials*ntap,nchan,nfoi,ntoi,cfg.precision);             end
            if csdflg, crsspctrm     = complex(zeros(ntrials*ntap,nchancmb,nfoi,ntoi,cfg.precision)); end
            if fftflg, fourierspctrm = complex(zeros(ntrials*ntap,nchan,nfoi,ntoi,cfg.precision));    end
            dimord    = 'rpttap_chan_freq_time';
          end
          
        case 'mtmfft'
          [spectrum,ntaper,foi] = specest_mtmfft(dat, time, options{:});
          nfoi = numel(foi);
          ntap = size(spectrum,1);
          if strcmp(cfg.calcdof,'yes')
            dof = zeros(ntrials,nfoi);
          end
          
          
          % allocate memory to output variables
          if keeprpt == 1 % cfg.keeptrials,'no' &&  cfg.keeptapers,'no'
            if powflg, powspctrm     = zeros(nchan,nfoi,cfg.precision);             end
            if csdflg, crsspctrm     = complex(zeros(nchancmb,nfoi,cfg.precision)); end
            if fftflg, fourierspctrm = complex(zeros(nchan,nfoi,cfg.precision));    end
            dimord    = 'chan_freq';
          elseif keeprpt == 2 % cfg.keeptrials,'yes' &&  cfg.keeptapers,'no'
            if powflg, powspctrm     = zeros(ntrials,nchan,nfoi,cfg.precision);             end
            if csdflg, crsspctrm     = complex(zeros(ntrials,nchancmb,nfoi,cfg.precision)); end
            if fftflg, fourierspctrm = complex(zeros(ntrials,nchan,nfoi,cfg.precision));    end
            dimord    = 'rpt_chan_freq';
          elseif keeprpt == 4 % cfg.keeptrials,'yes' &&  cfg.keeptapers,'yes'
            % FIXME this works only if all frequencies have the same number of tapers
            if powflg, powspctrm     = zeros(ntrials*ntap,nchan,nfoi,cfg.precision);             end
            if csdflg, crsspctrm     = complex(zeros(ntrials*ntap,nchancmb,nfoi,cfg.precision)); end
            if fftflg, fourierspctrm = complex(zeros(ntrials*ntap,nchan,nfoi,cfg.precision));    end
            dimord    = 'rpttap_chan_freq';
          end
        otherwise
          error('method %s is unknown', cfg.method);
      end % switch
    end % if itrial == 1
    
    
    
    
    
    % do the spectral decompisition of this trial and put it into output variables
    switch cfg.method
      
      case 'mtmconvol'
        if itrial ~= 1
          [spectrum,ntaper,foi,toi] = specest_mtmconvol(dat, time, 'timeoi', cfg.toi, options{:},'timwin',cfg.t_ftimwin);
          nfoi = numel(foi);
          ntoi = numel(toi);
          ntap = size(spectrum,1);
        end
        
        % get output in correct format
        % for now, there is a lot of redundancy, as each method has it's own case statement
        % when fully implemented, this can be cut down, perhaps in a separate switch, or perhaps as a time and a non-time if-loop
        foinumsmp = cfg.t_ftimwin .* data.fsample;
        foinumsmp = repmat(foinumsmp,[1, ntap, nchan, ntoi]);
        foinumsmp = permute(foinumsmp,[1 3 2 4]);
        if powflg
          powdum = 2.* abs(spectrum) .^ 2 ./ foinumsmp;
          if strcmp(cfg.taper, 'sine') % THIS IS NOT DONE IN THE MTMFFT CASE IN THE OLD IMPLEMENTATION, WHY?
            sinetapscale = zeros(ntap,nfoi);  % assumes fixed number of tapers
            for isinetap = 1:ntaper(1)  % assumes fixed number of tapers
              sinetapscale(isinetap,:) = (1 - (((isinetap - 1) ./ ntaper) .^ 2));
            end
            sinetapscale = reshape(repmat(sinetapscale,[1 1 nchan ntoi]),[ntap nchan nfoi ntoi]);
            powdum = powdum .* sinetapscale;
          end
          if keeprpt == 1 % cfg.keeptrials,'no' &&  cfg.keeptapers,'no'
            powspctrm = powspctrm + (reshape(nanmean(powdum,1),[nchan nfoi ntoi]) ./ ntrials);
          elseif keeprpt == 2  % cfg.keeptrials,'yes' &&  cfg.keeptapers,'no'
            powspctrm(itrial,:,:,:) = reshape(nanmean(powdum,1),[nchan nfoi ntoi]);
          elseif keeprpt == 4 % cfg.keeptrials,'yes' &&  cfg.keeptapers,'yes'
            rptind = reshape(1:ntrials * ntap,[ntap ntrials]);
            powspctrm(rptind(:,itrial),:,:,:) = reshape(powdum,[ntap nchan nfoi ntoi]);
          end
        end
        if fftflg
          fourierdum = spectrum .* sqrt(2 ./ foinumsmp); %cf Numercial Receipes 13.4.9
          if keeprpt == 1
            fourierspctrm = fourierspctrm + (reshape(nanmean(fourierdum,1),[nchan nfoi ntoi]) ./ ntrials);
          elseif keeprpt == 2
            fourierspctrm(itrial,:,:,:) = reshape(nanmean(fourierdum,1), [nchan nfoi ntoi]);
          elseif keeprpt == 4
            rptind = reshape(1:ntrials * ntap,[ntap ntrials]);
            fourierspctrm(rptind(:,itrial),:,:,:) = reshape(fourierdum,[ntap nchan nfoi ntoi]);
          end
        end
        if csdflg
          csddum = 2.* (spectrum(:,cutdatindcmb(:,1),:,:) .* conj(spectrum(:,cutdatindcmb(:,2),:,:))) ./ foinumsmp;
          if keeprpt == 1
            crsspctrm = crsspctrm + (reshape(nanmean(csddum,1),[nchancmb nfoi ntoi]) ./ ntrials);
          elseif keeprpt == 2
            crsspctrm(itrial,:,:,:) = reshape(nanmean(csddum,1), [nchancmb nfoi ntoi]);
          elseif keeprpt == 4
            rptind = reshape(1:ntrials * ntap,[ntap ntrials]);
            crsspctrm(rptind(:,itrial),:,:,:) = reshape(csddum,[ntap nchancmb nfoi ntoi]);
          end
        end
        
        if strcmp(cfg.calcdof,'yes')
          acttimboiind = ~isnan(squeeze(spectrum(1,1,:,:)));
          dof(itrial,:,:) = repmat(ntaper,[1, ntoi]) .* acttimboiind;
        end
        
        
        
      case 'mtmfft'
        if itrial ~= 1
          [spectrum,ntaper,foi] = specest_mtmfft(dat, time, options{:});
          nfoi = numel(foi);
          ntap = size(spectrum,1); % assumes fixed number of tapers
        end
        
        % get output in correct format
        % for now, there is a lot of redundancy, as each method has it's own case statement
        % when fully implemented, this can be cut down, perhaps in a separate switch, or perhaps as a time and a non-time if-loop
        foinumsmp = cfg.pad * data.fsample;
        foinumsmp = repmat(foinumsmp,[ntap, nchan, nfoi]);
        if powflg
          powdum = 2.* abs(spectrum) .^ 2 ./ foinumsmp;
          % In freqanalysis_mtmfft, sine tapers are NOT scaled (they are in mtmconvol). below code performs the scaling for cfg.method = mtmfft
          %           if strcmp(cfg.taper, 'sine')
          %             sinetapscale = zeros(ntap,nfoi);  % assumes fixed number of tapers
          %             for isinetap = 1:ntaper(1)  % assumes fixed number of tapers
          %               sinetapscale(isinetap,:) = (1 - (((isinetap - 1) ./ ntaper) .^ 2));
          %             end
          %             sinetapscale = reshape(repmat(sinetapscale,[1 1 nchan]),[ntap nchan nfoi]);
          %             powdum = powdum .* sinetapscale;
          %           end
          if keeprpt == 1 % cfg.keeptrials,'no' &&  cfg.keeptapers,'no'
            powspctrm = powspctrm + (reshape(nanmean(powdum,1),[nchan nfoi]) ./ ntrials);
          elseif keeprpt == 2  % cfg.keeptrials,'yes' &&  cfg.keeptapers,'no'
            powspctrm(itrial,:,:) = reshape(nanmean(powdum,1),[nchan nfoi]);
          elseif keeprpt == 4 % cfg.keeptrials,'yes' &&  cfg.keeptapers,'yes'
            rptind = reshape(1:ntrials * ntap,[ntap ntrials]);
            powspctrm(rptind(:,itrial),:,:) = reshape(powdum,[ntap nchan nfoi]);
          end
        end
        if fftflg
          fourierdum = spectrum .* sqrt(2 ./ foinumsmp); %cf Numercial Receipes 13.4.9
          if keeprpt == 1
            fourierspctrm = fourierspctrm + (reshape(nanmean(fourierdum,1),[nchan nfoi]) ./ ntrials);
          elseif keeprpt == 2
            fourierspctrm(itrial,:,:) = reshape(nanmean(fourierdum,1), [nchan nfoi]);
          elseif keeprpt == 4
            rptind = reshape(1:ntrials * ntap,[ntap ntrials]);
            fourierspctrm(rptind(:,itrial),:,:) = reshape(fourierdum,[ntap nchan nfoi]);
          end
        end
        if csdflg
          csddum = 2.* (spectrum(:,cutdatindcmb(:,1),:,:) .* conj(spectrum(:,cutdatindcmb(:,2),:,:))) ./ foinumsmp;
          if keeprpt == 1
            crsspctrm = crsspctrm + (reshape(nanmean(csddum,1),[nchancmb nfoi]) ./ ntrials);
          elseif keeprpt == 2
            crsspctrm(itrial,:,:) = reshape(nanmean(csddum,1), [nchancmb nfoi]);
          elseif keeprpt == 4
            rptind = reshape(1:ntrials * ntap,[ntap ntrials]);
            crsspctrm(rptind(:,itrial),:,:) = reshape(csddum,[ntap nchancmb nfoi]);
          end
        end
        
        if strcmp(cfg.calcdof,'yes')
          dof(itrial,:) = ntaper;
        end
        
        
        %       case 'wltconvol' % not testest yet
        %         [spectrum,foi,toi] = specest_wltconvol(dat, time, options{:});
        
        
      otherwise
        error('method %s is unknown', cfg.method);
    end % switch
    
    
    
    
    
    
  end % for ntrials
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% END: Main loop over trials
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  
  % set output variables
  freq = [];
  freq.label = data.label;
  freq.dimord = dimord;
  freq.freq   = foi;
  hasdc       = find(foi==0);
  if exist('toi','var')
    freq.time = toi;
  end
  if powflg
    % correct the 0 Hz bin if present, scaling with a factor of 2 is only appropriate for ~0 Hz
    if ~isempty(hasdc)
      if keeprpt>1
        powspctrm(:,:,hasdc,:) = powspctrm(:,:,hasdc,:)./2;      
      else
        powspctrm(:,hasdc,:) = powspctrm(:,hasdc,:)./2;
      end
    end
    freq.powspctrm = powspctrm;
  end
  if fftflg
    % correct the 0 Hz bin if present
    if ~isempty(hasdc)
      if keeprpt>1
        fourierspctrm(:,:,hasdc,:) = fourierspctrm(:,:,hasdc,:)./sqrt(2);      
      else
        fourierspctrm(:,hasdc,:) = fourierspctrm(:,hasdc,:)./sqrt(2);
      end
    end
    freq.fourierspctrm = fourierspctrm;
  end
  if csdflg
    % correct the 0 Hz bin if present
    if ~isempty(hasdc)
      if keeprpt>1
        crsspctrm(:,:,hasdc,:) = crsspctrm(:,:,hasdc,:)./2;      
      else
        crsspctrm(:,hasdc,:) = crsspctrm(:,hasdc,:)./2;
      end
    end
    freq.labelcmb  = cfg.channelcmb;
    freq.crsspctrm = crsspctrm;
  end
  if strcmp(cfg.calcdof,'yes');
    freq.dof = 2 .* dof;
  end;
  if strcmp(cfg.method,'mtmfft') && (keeprpt == 2 || keeprpt == 4)
    freq.cumsumcnt = trllength';
  end
  if keeprpt == 2,
    freq.cumtapcnt = repmat(ntaper(:)', [size(powspctrm,1) 1]);
  elseif keeprpt == 4,
    freq.cumtapcnt = repmat(ntaper(1), [size(fourierspctrm,1)./ntaper(1) 1]); % assumes fixed number of tapers
  end
  
  
  
  
  
  
  % accessing this field here is needed for the configuration tracking
  % by accessing it once, it will not be removed from the output cfg
  cfg.outputfile;
  
  % get the output cfg
  cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');
  
  try, freq.grad = data.grad; end   % remember the gradiometer array
  try, freq.elec = data.elec; end   % remember the electrode array
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






