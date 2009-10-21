function [freq] = freqanalysis_mtmfft(cfg, data);

% FREQANALYSIS_MTMFFT performs frequency analysis on any time series
% trial data using the 'multitaper method' (MTM) based on discrete
% prolate spheroidal sequences (Slepian sequences) as tapers. Alternatively,
% you can use conventional tapers (e.g. Hanning).
%
% Use as
%   [freq] = freqanalysis(cfg, data)
%
% The data should be organised in a structure as obtained from
% the PREPROCESSING function. The configuration should be according to
%   cfg.method     = method used for frequency or time-frequency decomposition
%                    see FREQANALYSIS for details
%   cfg.output     = 'pow'       return the power-spectra
%                    'powandcsd' return the power and the cross-spectra
%                    'fourier'   return the complex Fourier-spectra
%   cfg.taper      = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%
% For cfg.output='powandcsd', you should specify the channel combinations
% between which to compute the cross-spectra as cfg.channelcmb. Otherwise
% you should specify only the channels in cfg.channel.
%
%   cfg.channel    = Nx1 cell-array with selection of channels (default = 'all'),
%                    see CHANNELSELECTION for details
%   cfg.channelcmb = Mx2 cell-array with selection of channel pairs (default = {'all' 'all'}),
%                    see CHANNELCOMBINATION for details
%   cfg.foilim     = [begin end], frequency band of interest
%   cfg.tapsmofrq  = number, the amount of spectral smoothing through
%                    multi-tapering. Note that 4 Hz smoothing means
%                    plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
%   cfg.trials     = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.keeptrials = 'yes' or 'no', return individual trials or average (default = 'no')
%   cfg.keeptapers = 'yes' or 'no', return individual tapers or average (default = 'no')
%   cfg.pad        = number or 'maxperlen', length in seconds to which the data can be padded out (default = 'maxperlen')
%
% The padding will determine your spectral resolution. If you want to
% compare spectra from data pieces of different lengths, you should use
% the same cfg.pad for both, in order to spectrally interpolate them to
% the same spectral resolution.  Note that this will run very slow if you
% specify cfg.pad as maxperlen AND the number of samples turns out to have
% a large prime factor sum. This is because the FFTs will then be computed
% very inefficiently.
%
% See also FREQANALYSIS_MTMCONVOL, FREQANALYSIS_WLTCONVOL, FREQANALYSIS_TFR

% Undocumented local options
%   cfg.calcdof = 'yes'   calculate the degrees of freedom for every trial

% Copyright (c) 2003-2006, Pascal Fries, F.C. Donders Centre
%
% $Log: freqanalysis_mtmfft.m,v $
% Revision 1.45  2009/03/11 10:39:37  roboos
% more strict checking of cfg.pad
%
% Revision 1.44  2008/12/11 15:52:45  roboos
% fixed bug in normalisation of non-dpss tapers in case of non-rectangular data (i.e. data in which the length of the data is different over trials)
%
% Revision 1.43  2008/12/03 14:04:31  roboos
% only give warning about 1 taper for dpss
%
% Revision 1.42  2008/11/11 18:59:26  sashae
% added call to checkconfig at end of function (trackconfig and checksize)
%
% Revision 1.41  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.40  2008/05/13 15:36:09  roboos
% fixed potential bug in assessing the number of trials (when data.trial was column instead of row vector)a, now use numel instead of size
%
% Revision 1.39  2008/01/18 13:14:50  sashae
% added option for trial selection, updated documentation
%
% Revision 1.38  2007/09/23 14:10:09  erimar
% Replaced variable identifier as an argument to "exist" by its
% corresponding string.
%
% Revision 1.37  2007/08/06 15:02:03  roboos
% only changes in whitespace, no functional changes
%
% Revision 1.36  2007/03/27 11:00:28  ingnie
% deleted call to fixdimord because this is low level function
%
% Revision 1.35  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.34  2006/06/20 16:28:29  ingnie
% added consistent handling of cfg.channel and cfg.channelcmb
%
% Revision 1.33  2006/06/13 14:53:22  ingnie
% some change in white space defaults, added default cfg.channel = 'all'
%
% Revision 1.32  2006/06/06 16:57:51  ingnie
% updated documentation
%
% Revision 1.31  2006/05/04 14:27:00  roboos
% fixed a detail in the documentation
%
% Revision 1.30  2006/05/04 12:56:02  roboos
% updated documentation, removed powandfourier output
%
% Revision 1.29  2006/03/06 13:53:06  roboos
% pre-allocate sgncmbindx to save many memory allocation and copy operations
%
% Revision 1.28  2006/03/06 09:45:36  roboos
% fixed the callback detection for octave, added sine tapet (thanks to Tom)
% changed some | into ||
%
% Revision 1.27  2006/02/28 12:23:10  erimar
% Added configuratio option cfg.calcdof to request for the calculation of the degrees of freedom.
%
% Revision 1.26  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.25  2006/02/14 08:41:41  roboos
% replaced strmatch with find(strcmp) to speed up the handling of channelcombinations in Octave
%
% Revision 1.24  2006/02/07 20:08:22  roboos
% changed all occurences of a dimord with chancmb (was previous sgncmb) into chan
%
% Revision 1.23  2006/02/01 12:26:00  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.22  2005/09/22 14:29:20  jansch
% correct normalisation for tapers, different than dpss
%
% Revision 1.21  2005/08/23 12:28:36  jansch
% be sure that keeptapers and keeptrials = 'yes' for fourier-output
%
% Revision 1.20  2005/08/20 08:13:47  roboos
% do not try to normalise powspctrm if output=foerier
%
% Revision 1.19  2005/08/19 08:11:45  roboos
% implemented local subfunction that ensures that the input arguments for dpss are double precision
%
% Revision 1.18  2005/08/19 08:10:08  roboos
% changed default output: powandcsd if channelcmb specified, pow otherwise
%
% Revision 1.17  2005/08/16 07:17:11  jansch
% Implemented the undocumented option of outputting the fourierspectrum only.
% Eventually this is intended to replace the 'powandfourier'-option altogether.
%
% Revision 1.16  2005/08/15 13:30:11  jansch
% changed the normalisation of the single-trial fourierspectra so that it is
% consistent with how it is done for the powerspectra and the cross-spectra.
%
% Revision 1.15  2005/08/09 10:15:15  jansch
% implemented pre-allocation of memory in the case of keeptapers. changed the
% cumsumcnt and cumtapcnt in the case of keeptapers. TAKE CARE, this is not
% necessarily backwardcompatible. cumsumcnt and cumtapcnt are now expressed
% in terms of the original number of trials, and not anymore on the total
% number of tapers.
%
% Revision 1.14  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.13  2005/06/09 07:12:14  jansch
% fixed bug in assigning cumtapcnt to output in the case of keeptaper
%
% Revision 1.12  2005/06/08 09:48:28  jansch
% added the tapercounter for the keeptapers-option
%
% Revision 1.11  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.10  2005/05/10 07:39:33  jansch
% fixed small bug of assignment of number of tapers per trial in the case of
% equal-length trials
%
% Revision 1.9  2005/02/16 08:43:02  jansch
% changed f_t into fourier
%
% Revision 1.8  2005/02/15 11:02:25  jansch
% added the (hidden) option to keep the complex fourier-transform in the output.
% therefore, cfg.output has to be set to 'powandfourier'. this yields, in addition to
% the well-known powspctrm, an fourierspctrm, which can be used to extract the
% phase from the individual trials by using the angle command.
%
% Revision 1.7  2005/01/19 08:42:42  jansch
% removed obsolete code for generating and checking channelcombinations
% cleaned up handling of sgnindx/sgncmbindx
% ensured that power is computed for all channels that are in channelcmb
%
% Revision 1.6  2005/01/18 15:11:39  roboos
% Cleaned up configuration for sgn/sgncmb, now exclusively using channel/channelcmb which is consistent with rest of fieldtrip and freqanalysis documentation. Also the output now only contains freq.label/labelcmb and not any more the old sgn/sgncmb.
%
% Revision 1.5  2005/01/17 14:51:16  roboos
% implemented tapering with single window according to user specification
% (cfg.taper) to allow hanning and other windows to be used easily.
% The default is to multitaper with a dpss sequence.
%
% Revision 1.4  2004/11/18 16:51:58  roboos
% updated the help: pointed to CHANNELCOMBINATION for cfg.channelcmb
%
% Revision 1.3  2004/10/28 09:00:33  roboos
% fixed bug in caller function detection (for matlab 6.1 and 6.5)
%
% Revision 1.2  2004/10/28 07:21:46  roboos
% added check to ensure that the FREQANALYSIS wrapper is started instead of the
% FERQANALYSIS_xxx subfunctions
%
% Revision 1.1  2004/09/21 12:02:04  marsie
% the mtmfft method of multitaperanalysis.m has been moved to this seperate function
%
% Revision 1.3  2004/09/16 11:40:34  marsie
% fixed bug in 'mtmfft' that slowed computation down w/o causing bad results
%
% Revision 1.2  2004/08/31 14:36:30  marsie
% added support for data segments, which are to short to place a taper
%
% Revision 1.1  2004/08/25 17:38:28  marsie
% initial CVS commit
% this version is based on the original freqanalysis including the following changes
% and improvements:
% - the DFT method is no longer supported
% - methods are renamed with the prefix 'mtm': 'mtmfft' and 'mtmconvol'
% - the last taper returned by 'dpss' is not used
% - the function checks for the consistency of data length and spectral smoothing
% - improved performance by vectorization of loops
% - output of 'mtmfft' is correctly dimensioned for single frequency output
%

fieldtripdefs

% ensure that this function is started as a subfunction of the FREQANALYSIS wrapper
if ~exist('OCTAVE_VERSION')
  [s, i] = dbstack;
  if length(s)>1
    [caller_path, caller_name, caller_ext] = fileparts(s(2).name);
  else
    caller_path = '';
    caller_name = '';
    caller_ext  = '';
  end
  % evalin('caller', 'mfilename') does not work for Matlab 6.1 and 6.5
  if ~strcmp(caller_name, 'freqanalysis')
    error(['you should call FREQANALYSIS, instead of the ' upper(mfilename) ' subfunction']);
  end
end

% set all the defaults
if ~isfield(cfg, 'method'),        cfg.method     = 'mtmfft';               end
if ~isfield(cfg, 'keeptapers'),    cfg.keeptapers = 'no';                   end
if ~isfield(cfg, 'keeptrials'),    cfg.keeptrials = 'no';                   end
if ~isfield(cfg, 'calcdof'),       cfg.calcdof    = 'no';                   end
if ~isfield(cfg, 'pad'),           cfg.pad        = 'maxperlen';            end
if ~isfield(cfg, 'taper'),         cfg.taper      = 'dpss';                 end
if ~isfield(cfg, 'channel'),       cfg.channel    = 'all';                  end
if ~isfield(cfg, 'foilim'),        cfg.foilim     = [0 data.fsample/2];     end
if ~isfield(cfg, 'output'),
  if isfield(cfg, 'channelcmb') && ~isempty(cfg.channelcmb)
    cfg.output     = 'powandcsd';
  else
    cfg.output     = 'pow';
  end
end

if strcmp(cfg.output, 'fourier'),
  cfg.keeptrials = 'yes';
  cfg.keeptapers = 'yes';
end

if ~strcmp(cfg.method,'mtmfft')
  error('unsupported method');
end

% setting a flag (csdflg) that determines whether this routine outputs
% only power-spectra or power-spectra and cross-spectra?
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
  error('unsupported value for cfg.method');
end

if ~isfield(cfg, 'channelcmb') && csdflg
  %set the default for the channelcombination
  cfg.channelcmb = {'all' 'all'};
elseif isfield(cfg, 'channelcmb') && ~csdflg
  % no cross-spectrum needs to be computed, hence remove the combinations from cfg
  cfg = rmfield(cfg, 'channelcmb');
end

% ensure that channelselection and selection of channelcombinations is
% perfomed consistently
cfg.channel = channelselection(cfg.channel, data.label);
if isfield(cfg, 'channelcmb')
  cfg.channelcmb = channelcombination(cfg.channelcmb, data.label);
end

% determine the corresponding indices of all channels
sgnindx     = match_str(data.label, cfg.channel);
numsgn      = size(sgnindx,1);
if csdflg
  % determine the corresponding indices of all channel combinations
  sgncmbindx = zeros(size(cfg.channelcmb));
  for k=1:size(cfg.channelcmb,1)
    sgncmbindx(k,1) = find(strcmp(cfg.channelcmb(k,1), data.label));
    sgncmbindx(k,2) = find(strcmp(cfg.channelcmb(k,2), data.label));
    % this works the same, but is much slower in Octave
    % sgncmbindx(k,1) = strmatch(cfg.channelcmb(k,1), data.label, 'exact');
    % sgncmbindx(k,2) = strmatch(cfg.channelcmb(k,2), data.label, 'exact');
  end

  numsgncmb   = size(sgncmbindx,1);
  sgnindx     = unique([sgnindx(:); sgncmbindx(:)]);
  numsgn      = length(sgnindx);

  cutdatindcmb = zeros(size(sgncmbindx));
  for sgnlop = 1:numsgn
    cutdatindcmb(find(sgncmbindx == sgnindx(sgnlop))) = sgnlop;
  end
end

% if rectan is 1 it means that trials are of equal lengths
numper = numel(data.trial);
rectan = 1;
for perlop = 1:numper
  numdatbnsarr(perlop,1) = size(data.trial{perlop},2);
end
rectan = all(numdatbnsarr==numdatbnsarr(1));

% if cfg.pad is 'maxperlen', this is realized here:
if isequal(cfg.pad, 'maxperlen')
  cfg.pad = max(numdatbnsarr,[],1) ./ data.fsample;
else
  % check that the specified padding is not too short
  if cfg.pad<(max(numdatbnsarr,[],1)/data.fsample)
    error('the padding that you specified is shorter than the longest trial in the data');
  end
end
numsmp = ceil(cfg.pad .* data.fsample); % this used to be "cfg.pad .* data.fsample"

% keeping trials and/or tapers?
if strcmp(cfg.keeptrials,'no') &&  strcmp(cfg.keeptapers,'no')
  keep = 1;
elseif strcmp(cfg.keeptrials,'yes') &&  strcmp(cfg.keeptapers,'no')
  keep = 2;
elseif strcmp(cfg.keeptrials,'no') &&  strcmp(cfg.keeptapers,'yes')
  error('There is no support for keeping tapers WITHOUT KEEPING TRIALS.');
elseif strcmp(cfg.keeptrials,'yes') &&  strcmp(cfg.keeptapers,'yes')
  keep = 4;
end

% calculating degrees of freedom
calcdof = strcmp(cfg.calcdof,'yes');

% doing the computation
boilim  = round(cfg.foilim ./ (data.fsample ./ numsmp)) + 1;
boi     = boilim(1):boilim(2);
numboi  = size(boi,2);
foi     = (boi-1) ./ cfg.pad;

if keep == 1
  if powflg, powspctrm     = zeros(numsgn,numboi);              end
  if csdflg, crsspctrm     = complex(zeros(numsgncmb,numboi));  end
  if fftflg, fourierspctrm = complex(zeros(numsgn,numboi));     end
  dimord    = 'chan_freq';
elseif keep == 2
  if powflg, powspctrm     = zeros(numper,numsgn,numboi);             end
  if csdflg, crsspctrm     = complex(zeros(numper,numsgncmb,numboi)); end
  if fftflg, fourierspctrm = complex(zeros(numper,numsgn,numboi));    end
  dimord    = 'rpt_chan_freq';
elseif keep == 4
  if rectan == 1, % compute the amount of memory needed to collect the results
    numdatbns = numdatbnsarr(1,1);
    if strcmp(cfg.taper, 'dpss'),
      % ensure that the input arguments are double precision
      tap = double_dpss(numdatbns, numdatbns*(cfg.tapsmofrq./data.fsample))';
    elseif strcmp(cfg.taper, 'sine')
      tap = sine_taper(numdatbns, numdatbns*(cfg.tapsmofrq./data.fsample))';
    else
      tap(2,:) = nan;
    end
    numtap = size(tap,1)-1;
    numrpt = numtap.*numper;
  elseif rectan == 0,
    numrpt = 0;
    for perlop = 1:numper
      numdatbns = numdatbnsarr(perlop,1);
      if strcmp(cfg.taper, 'dpss'),
        % ensure that the input arguments are double precision
        tap = double_dpss(numdatbns, numdatbns*(cfg.tapsmofrq./data.fsample))';
      elseif strcmp(cfg.taper, 'sine')
        tap = sine_taper(numdatbns, numdatbns*(cfg.tapsmofrq./data.fsample))';
      else
        tap(2,:) = nan;
      end
      numtap = size(tap,1)-1;
      numrpt = numrpt + numtap;
    end
  end
  if powflg, powspctrm     = zeros(numrpt,numsgn,numboi);             end
  if csdflg, crsspctrm     = complex(zeros(numrpt,numsgncmb,numboi)); end
  if fftflg, fourierspctrm = complex(zeros(numrpt,numsgn,numboi));    end
  cnt = 0;
  dimord    = 'rpttap_chan_freq';
end

% these count the number of tapers
cumsumcnt = zeros(numper,1);
cumtapcnt = zeros(numper,1);

if rectan == 1
  % trials are of equal length, compute the set of tapers only once
  numdatbns = numdatbnsarr(1,1);
  if strcmp(cfg.taper, 'dpss')
    % create a sequence of DPSS (Slepian) tapers
    % ensure that the input arguments are double precision
    tap = double_dpss(numdatbns,numdatbns*(cfg.tapsmofrq./data.fsample))';
  elseif strcmp(cfg.taper, 'sine')
    tap = sine_taper(numdatbns, numdatbns*(cfg.tapsmofrq./data.fsample))';
  else
    % create a single taper according to the window specification as a
    % replacement for the DPSS (Slepian) sequence
    tap = window(cfg.taper, numdatbns)';
    tap = tap./norm(tap);
    % freqanalysis_mtmfft always throws away the last taper of the Slepian sequence, so add a dummy taper
    tap(2,:) = nan;
  end
  numtap = size(tap,1) - 1;
  if keep == 2 || calcdof
    cumtapcnt(:) = numtap;
  end
  if (numtap < 1)
    error(sprintf(...
      'datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',...
      numdatbns/data.fsample, cfg.tapsmofrq, data.fsample/numdatbns));
  elseif (numtap < 2) && strcmp(cfg.taper, 'dpss')
    fprintf('WARNING: using only one taper for specified smoothing\n');
  end
  pad = zeros(1,numsmp - numdatbns);
end

for perlop = 1:numper
  fprintf('processing trial %d, ', perlop);
  if keep == 2
    cnt = perlop;
    cumsumcnt(cnt,1) = numdatbnsarr(perlop,1);
  end
  if rectan == 0
    % trials are not of equal length, compute the set of tapers for this trial
    numdatbns = numdatbnsarr(perlop,1);
    if strcmp(cfg.taper, 'dpss')
      % create a sequence of DPSS (Slepian) tapers
      % ensure that the input arguments are double precision
      tap = double_dpss(numdatbns,numdatbns*(cfg.tapsmofrq./data.fsample))';
    elseif strcmp(cfg.taper, 'sine')
      tap = sine_taper(numdatbns, numdatbns*(cfg.tapsmofrq./data.fsample))';
    else
      % create a single taper according to the window specification as a
      % replacement for the DPSS (Slepian) sequence
      tap = window(cfg.taper, numdatbns)';
      tap = tap./norm(tap);
      % freqanalysis_mtmfft always throws away the last taper of the Slepian sequence, so add a dummy taper
      tap(2,:) = nan;
    end
    numtap = size(tap,1) - 1;
    if keep == 2 || calcdof
      cnt = perlop;
      cumtapcnt(cnt,1) = numtap;
    end
    if (numtap < 1)
      error(sprintf(...
        'datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',...
        numdatbns/data.fsample, cfg.tapsmofrq, data.fsample/numdatbns));
    elseif (numtap < 2) && strcmp(cfg.taper, 'dpss')
      fprintf('WARNING: using only one taper for specified smoothing\n');
    end
    pad = zeros(1,numsmp - numdatbns);
  end
  for taplop = 1:numtap
    if keep == 4
      cnt = cnt+1;
      cumsumcnt(perlop,1) = numdatbnsarr(perlop,1);
      cumtapcnt(perlop,1) = numtap;
    end
    if calcdof
      cumtapcnt(perlop,1) = numtap;
    end

    autspctrmacttap = complex(zeros(numsgn,numboi));
    for sgnlop = 1:numsgn
      dum = fft([data.trial{perlop}(sgnindx(sgnlop),:) .* tap(taplop,:) , ...
        pad],[],2);
      autspctrmacttap(sgnlop,:) = dum(boi);
    end
    if taplop == 1
      fprintf('nfft: %d samples, taper length: %d samples, %d tapers\n',length(dum),size(tap,2),numtap);
    end
    if powflg
      powdum = 2 .* (autspctrmacttap .* conj(autspctrmacttap)) ./ numsmp; %cf Numercial Receipes 13.4.9
      if keep == 1
        powspctrm(:,:) = powspctrm(:,:) + (powdum ./ numtap);
      elseif keep == 2
        powspctrm(cnt,:,:) = powspctrm(cnt,:,:) + (permute(powdum,[3,1,2]) ./ numtap);
      elseif keep == 4
        powspctrm(cnt,:,:) = powdum;
      end
    end
    if fftflg
      fourierdum = (autspctrmacttap) .* sqrt(2 ./ numsmp); %cf Numercial Receipes 13.4.9
      if keep == 1
        fourierspctrm(:,:) = fourierspctrm(:,:) + (fourierdum ./ numtap);
      elseif keep == 2
        fourierspctrm(cnt,:,:) = fourierspctrm(cnt,:,:) + (permute(fourierdum,[3,1,2]) ./ numtap);
      elseif keep == 4
        fourierspctrm(cnt,:,:) = fourierdum;
      end
    end
    if csdflg
      csddum = 2.* (autspctrmacttap(cutdatindcmb(:,1),:) .* ...
        conj(autspctrmacttap(cutdatindcmb(:,2),:))) ./ numsmp;
      if keep == 1
        crsspctrm(:,:) = crsspctrm(:,:) + csddum ./ numtap;
      elseif keep == 2
        crsspctrm(cnt,:,:) = crsspctrm(cnt,:,:) + permute(csddum,[3,1,2]) ./ numtap;
      elseif keep == 4
        crsspctrm(cnt,:,:) = csddum;
      end
    end
  end % taplop
end % perlop
if keep ==1
  if powflg, powspctrm = powspctrm ./ numper; end
  if csdflg, crsspctrm = crsspctrm ./ numper; end
end

% collect the results
freq.label      = data.label(sgnindx);
freq.dimord     = dimord;
freq.freq       = foi;
if powflg
  freq.powspctrm  = powspctrm;
end
if fftflg
  freq.fourierspctrm  = fourierspctrm;
end

if csdflg
  freq.labelcmb   = cfg.channelcmb;
  freq.crsspctrm  = crsspctrm;
end

if strcmp(cfg.method,'mtmfft') && (keep == 2 || keep == 4)
  freq.cumsumcnt = cumsumcnt;
end

if strcmp(cfg.method,'mtmfft') && (keep == 2 || keep == 4)
  freq.cumtapcnt = cumtapcnt;
end

if calcdof
  freq.dof=2*repmat(cumtapcnt,[1,numboi]);
end;

try, freq.grad = data.grad; end   % remember the gradiometer array
try, freq.elec = data.elec; end   % remember the electrode array

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
cfg.version.id = '$Id: freqanalysis_mtmfft.m,v 1.45 2009/03/11 10:39:37 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
freq.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION ensure that the first two input arguments are of double
% precision this prevents an instability (bug) in the computation of the
% tapers for Matlab 6.5 and 7.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tap] = double_dpss(a, b, varargin);
tap = dpss(double(a), double(b), varargin{:});

