function [freq] = freqanalysis_mtmconvol(cfg, data);

% FREQANALYSIS_MTMCONVOL performs time-frequency analysis on any time series trial data
% using the 'multitaper method' (MTM) based on Slepian sequences as tapers. Alternatively,
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
%   cfg.foi        = vector 1 x numfoi, frequencies of interest
%   cfg.t_ftimwin  = vector 1 x numfoi, length of time window (in seconds)
%   cfg.tapsmofrq  = vector 1 x numfoi, the amount of spectral smoothing through
%                    multi-tapering. Note that 4 Hz smoothing means
%                    plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
%   cfg.toi        = vector 1 x numtoi, the times on which the analysis windows
%                    should be centered (in seconds)
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
% An asymmetric taper usefull for TMS is used when cfg.taper='alpha' and corresponds to
%   W. Kyle Mitchell, Mark R. Baker & Stuart N.  Baker. Muscle Responses to
%   Transcranial Stimulation Depend on Background Oscillatory Activity.
%   published online Jul 12, 2007 J. Physiol. 
%
% See also FREQANALYSIS

% undocumented experimental options
%   cfg.calcdof = 'yes'   calculate the degrees of freedom for every trial

% Copyright (c) 2003,2004-2006 F.C. Donders Centre
%
% $Log: freqanalysis_mtmconvol.m,v $
% Revision 1.44  2009/03/11 10:39:37  roboos
% more strict checking of cfg.pad
%
% Revision 1.43  2008/11/11 18:59:26  sashae
% added call to checkconfig at end of function (trackconfig and checksize)
%
% Revision 1.42  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.41  2008/05/13 13:53:22  roboos
% added some documentation on alpha taper, fixed potential bug in assessing number of trials
%
% Revision 1.40  2008/02/28 09:44:30  roboos
% round the number of padded samples to an integer
%
% Revision 1.39  2008/01/18 13:14:50  sashae
% added option for trial selection, updated documentation
%
% Revision 1.38  2007/08/06 15:06:00  roboos
% implemented support for asymmetric alpha tapers
% changed the pre-allocation of some small arrays
% some changes in whitespace and removed some "..." to break up lines
%
% Revision 1.37  2007/03/27 11:00:28  ingnie
% deleted call to fixdimord because this is low level function
%
% Revision 1.36  2007/02/19 09:33:29  jansch
% fixed inappropriate behavour of freq.cumtapcnt in the case of fourierspectra as output
%
% Revision 1.35  2007/02/13 14:17:51  roboos
% made zero padding more memory efficient
%
% Revision 1.34  2006/11/30 07:51:12  jansch
% changed cumtapcnt into dof throughout the code. added cumtapcnt to output
%
% Revision 1.33  2006/10/19 15:32:55  roboos
% updated documentation
%
% Revision 1.32  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.31  2006/06/20 16:27:52  ingnie
% updated documentation, added consistent handling of cfg.channel and cfg.channelcmb
%
% Revision 1.30  2006/06/13 14:53:21  ingnie
% some change in white space defaults, added default cfg.channel = 'all'
%
% Revision 1.29  2006/06/06 16:57:51  ingnie
% updated documentation
%
% Revision 1.28  2006/05/22 14:47:04  jansch
% fixed a bug in check for different amount of tapers for different frequencies
%
% Revision 1.27  2006/05/03 09:52:22  jansch
% updated error-message when checking for keeptrials and keeptapers in combination
% with fourier as an output
%
% Revision 1.26  2006/03/06 13:53:06  roboos
% pre-allocate sgncmbindx to save many memory allocation and copy operations
%
% Revision 1.25  2006/03/06 09:45:36  roboos
% fixed the callback detection for octave, added sine tapet (thanks to Tom)
% changed some | into ||
%
% Revision 1.24  2006/02/28 12:25:09  erimar
% Added the configuration option cfg.calcdof to request for the calculation of the degrees of freedom.
%
% Revision 1.23  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.22  2006/02/07 20:08:22  roboos
% changed all occurences of a dimord with chancmb (was previous sgncmb) into chan
%
% Revision 1.21  2006/02/01 12:26:00  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.20  2006/01/30 14:04:21  jansch
% included the option to keep fourier-data. it only works when the number of
% tapers per frequency are equal across frequencies
%
% Revision 1.19  2005/10/13 11:33:36  roboos
% fixed bug in computation of cfg.pad in case of 'maxperlength' (thanks to Floris)
%
% Revision 1.18  2005/09/22 14:29:20  jansch
% correct normalisation for tapers, different than dpss
%
% Revision 1.17  2005/08/23 12:28:36  jansch
% be sure that keeptapers and keeptrials = 'yes' for fourier-output
%
% Revision 1.16  2005/08/19 08:11:45  roboos
% implemented local subfunction that ensures that the input arguments for dpss are double precision
%
% Revision 1.15  2005/08/15 13:30:11  jansch
% changed the normalisation of the single-trial fourierspectra so that it is
% consistent with how it is done for the powerspectra and the cross-spectra.
%
% Revision 1.14  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.13  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.12  2005/05/04 07:31:36  roboos
% fixed bug in non-dpss taper, had to be transposed
%
% Revision 1.11  2005/02/16 09:21:24  jansch
% inserted the cfg.output-option 'powandfourier'. in addition to outputting
% the powerspectrum, the complex discrete fourier spectrum is returned,
% allowing for computation of phases.
%
% Revision 1.10  2005/01/19 08:42:42  jansch
% removed obsolete code for generating and checking channelcombinations
% cleaned up handling of sgnindx/sgncmbindx
% ensured that power is computed for all channels that are in channelcmb
%
% Revision 1.9  2005/01/18 15:11:38  roboos
% Cleaned up configuration for sgn/sgncmb, now exclusively using channel/channelcmb which is consistent with rest of fieldtrip and freqanalysis documentation. Also the output now only contains freq.label/labelcmb and not any more the old sgn/sgncmb.
%
% Revision 1.8  2005/01/17 14:51:16  roboos
% implemented tapering with single window according to user specification
% (cfg.taper) to allow hanning and other windows to be used easily.
% The default is to multitaper with a dpss sequence.
%
% Revision 1.7  2004/12/20 14:52:56  roboos
% changed rounding off of timboi
%
% Revision 1.6  2004/12/20 14:47:03  roboos
% fixed rounding off error in sample number timboi
%
% Revision 1.5  2004/12/20 13:03:14  jansch
% fixed bug that incorrectly handled missing data as zeros
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
% Revision 1.1  2004/09/21 12:04:02  marsie
% the mtmconvol method of multitaperanalysis.m has been moved to this seperate function
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
if ~isfield(cfg, 'method'),        cfg.method     = 'mtmconvol';  end
if ~isfield(cfg, 'keeptapers'),    cfg.keeptapers = 'no';         end
if ~isfield(cfg, 'keeptrials'),    cfg.keeptrials = 'no';         end
if ~isfield(cfg, 'calcdof'),       cfg.calcdof    = 'no';         end
if ~isfield(cfg, 'output'),        cfg.output     = 'powandcsd';  end
if ~isfield(cfg, 'pad'),           cfg.pad        = 'maxperlen';  end
if ~isfield(cfg, 'taper'),         cfg.taper      = 'dpss';       end
if ~isfield(cfg, 'channel'),       cfg.channel    = 'all';        end
if strcmp(cfg.output, 'fourier'),
  cfg.keeptrials = 'yes';
  cfg.keeptapers = 'yes';
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
  error('Unrecognized output required');
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
    sgncmbindx(k,1) = strmatch(cfg.channelcmb(k,1), data.label, 'exact');
    sgncmbindx(k,2) = strmatch(cfg.channelcmb(k,2), data.label, 'exact');
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
numper       = numel(data.trial);
numdatbnsarr = zeros(numper, 1);
for perlop = 1:numper
  numdatbnsarr(perlop) = size(data.trial{perlop},2);
end
rectan = all(numdatbnsarr==numdatbnsarr(1));

% if cfg.pad is 'maxperlen', this is realized here:
% first establish where the first possible sample is
min_smp = min(data.offset);
% then establish where the last possible sample is
max_smp = max(numdatbnsarr(:)+data.offset(:));
if isequal(cfg.pad, 'maxperlen')
  % pad the data from the first possible to last possible sample
  cfg.pad = (max_smp-min_smp) ./ data.fsample;
else
  % check that the specified padding is not too short
  if cfg.pad<((max_smp-min_smp)/data.fsample)
    error('the padding that you specified is shorter than the longest trial in the data');
  end
end
  clear min_smp max_smp
numsmp = round(cfg.pad .* data.fsample);

% keeping trials and/or tapers?
if strcmp(cfg.keeptrials,'no') &&  strcmp(cfg.keeptapers,'no')
  keep = 1;
elseif strcmp(cfg.keeptrials,'yes') &&  strcmp(cfg.keeptapers,'no')
  keep = 2;
elseif strcmp(cfg.keeptrials,'no') &&  strcmp(cfg.keeptapers,'yes')
  error('There is currently no support for keeping tapers WITHOUT KEEPING TRIALS.');
elseif strcmp(cfg.keeptrials,'yes') &&  strcmp(cfg.keeptapers,'yes')
  keep = 4;
end

if strcmp(cfg.keeptrials,'yes') && strcmp(cfg.keeptapers,'yes')
  if ~strcmp(cfg.output, 'fourier'),
    error('Keeping trials AND tapers is only possible with fourier as the output.');
  elseif strcmp(cfg.taper, 'dpss') && ~(all(cfg.tapsmofrq==cfg.tapsmofrq(1)) && all(cfg.t_ftimwin==cfg.t_ftimwin(1))),
    error('Currently you can only keep trials AND tapers, when using the number of tapers per frequency is equal across frequency');
  end
end

if strcmp(cfg.taper, 'alpha') && ~all(cfg.t_ftimwin==cfg.t_ftimwin(1))
  error('you can only use alpha tapers with an cfg.t_ftimwin that is equal for all frequencies');
end

minoffset = min(data.offset);
timboi = round(cfg.toi .* data.fsample - minoffset);
toi    = round(cfg.toi .* data.fsample) ./ data.fsample;
numtoi = length(cfg.toi);
numfoi = length(cfg.foi);
numtap = zeros(numfoi,1);

% calculating degrees of freedom
calcdof = strcmp(cfg.calcdof,'yes');
if calcdof
  dof = zeros(numper,numfoi,numtoi);
end;

% compute the tapers and their fft
knlspctrmstr = cell(numfoi,1);
for foilop = 1:numfoi
  acttapnumsmp = round(cfg.t_ftimwin(foilop) .* data.fsample);
  if strcmp(cfg.taper, 'dpss')
    % create a sequence of DPSS (Slepian) tapers, ensure that the input arguments are double
    tap = double_dpss(acttapnumsmp, acttapnumsmp .* (cfg.tapsmofrq(foilop)./data.fsample));
  elseif strcmp(cfg.taper, 'sine')
    tap = sine_taper(acttapnumsmp, acttapnumsmp .* (cfg.tapsmofrq(foilop)./data.fsample));
  elseif strcmp(cfg.taper, 'alpha')
    tap = alpha_taper(acttapnumsmp, cfg.foi(foilop)./data.fsample);
    tap = tap./norm(tap);
    % freqanalysis_mtmconvol always throws away the last taper of the Slepian sequence, so add a dummy taper
    tap(:,2) = nan;
  else
    % create a single taper according to the window specification as a replacement for the DPSS (Slepian) sequence
    tap = window(cfg.taper, acttapnumsmp);
    tap = tap./norm(tap);
    % freqanalysis_mtmconvol always throws away the last taper of the Slepian sequence, so add a dummy taper
    tap(:,2) = nan;
  end
  numtap(foilop) = size(tap,2)-1;
  if (numtap(foilop) < 1)
    error(sprintf('%.3f Hz : datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz', cfg.foi(foilop), acttapnumsmp/data.fsample, cfg.tapsmofrq(foilop), data.fsample/acttapnumsmp));
  elseif (numtap(foilop) < 2) && strcmp(cfg.taper, 'dpss')
    fprintf('%.3f Hz : WARNING - using only one taper for specified smoothing\n',cfg.foi(foilop));
  end
  ins = ceil(numsmp./2) - floor(acttapnumsmp./2);
  prezer = zeros(ins,1);
  pstzer = zeros(numsmp - ((ins-1) + acttapnumsmp)-1,1);
  ind    = (0:acttapnumsmp-1)' .* ((2.*pi./data.fsample) .* cfg.foi(foilop));
  knlspctrmstr{foilop} = complex(zeros(numtap(foilop),numsmp));
  for taplop = 1:numtap(foilop)
    try
      % construct the complex wavelet
      coswav  = vertcat(prezer,tap(:,taplop).*cos(ind),pstzer);
      sinwav  = vertcat(prezer,tap(:,taplop).*sin(ind),pstzer);
      wavelet = complex(coswav, sinwav);
      % store the fft of the complex wavelet
      knlspctrmstr{foilop}(taplop,:) = fft(wavelet,[],1)';
      global fb
      if ~isempty(fb) && fb
        % plot the wavelet for debugging
        figure
        plot(tap(:,taplop).*cos(ind), 'r'); hold on
        plot(tap(:,taplop).*sin(ind), 'g');
        plot(tap(:,taplop)          , 'b');
        title(sprintf('taper %d @ %g Hz', taplop, cfg.foi(foilop)));
        drawnow
      end
    end
  end
end

if keep == 1
  if powflg, powspctrm     = zeros(numsgn,numfoi,numtoi);             end
  if csdflg, crsspctrm     = complex(zeros(numsgncmb,numfoi,numtoi)); end
  if fftflg, fourierspctrm = complex(zeros(numsgn,numfoi,numtoi));    end
  cntpertoi = zeros(numfoi,numtoi);
  dimord    = 'chan_freq_time';
elseif keep == 2
  if powflg, powspctrm     = zeros(numper,numsgn,numfoi,numtoi);             end
  if csdflg, crsspctrm     = complex(zeros(numper,numsgncmb,numfoi,numtoi)); end
  if fftflg, fourierspctrm = complex(zeros(numper,numsgn,numfoi,numtoi));    end
  dimord    = 'rpt_chan_freq_time';
elseif keep == 4
  % FIXME this works only if all frequencies have the same number of tapers
  if powflg, powspctrm     = zeros(numper*numtap(1),numsgn,numfoi,numtoi);             end
  if csdflg, crsspctrm     = complex(zeros(numper*numtap(1),numsgncmb,numfoi,numtoi)); end
  if fftflg, fourierspctrm = complex(zeros(numper*numtap(1),numsgn,numfoi,numtoi));    end
  cnt = 0;
  dimord    = 'rpttap_chan_freq_time';
end

for perlop = 1:numper
  fprintf('processing trial %d: %d samples\n', perlop, numdatbnsarr(perlop,1));
  if keep == 2
    cnt = perlop;
  end
  numdatbns = numdatbnsarr(perlop,1);
  % prepad = zeros(1,data.offset(perlop) - minoffset);
  % pstpad = zeros(1,minoffset + numsmp - (data.offset(perlop) + numdatbns));
  % datspctra = complex(zeros(numsgn,numsmp));
  % for sgnlop = 1:numsgn
  %   datspctra(sgnlop,:) = fft([prepad, data.trial{perlop}(sgnindx(sgnlop),:), ...
  %       pstpad],[],2);
  % end
  prepad = zeros(numsgn,data.offset(perlop) - minoffset);
  pstpad = zeros(numsgn,minoffset + numsmp - (data.offset(perlop) + numdatbns));
  tmp = data.trial{perlop}(sgnindx,:);
  tmp = [prepad tmp pstpad];
  % avoid the use of a 3rd input argument to facilitate compatibility with star-P
  % use explicit transpose, to avoid complex conjugate transpose
  datspctra = transpose(fft(transpose(tmp)));
  for foilop = 1:numfoi
    fprintf('processing frequency %d (%.2f Hz), %d tapers\n', foilop,cfg.foi(foilop),numtap(foilop));
    actfoinumsmp    = cfg.t_ftimwin(foilop) .* data.fsample;
    acttimboiind    = find(timboi >= (-minoffset + data.offset(perlop) + (actfoinumsmp ./ 2)) & timboi <  (-minoffset + data.offset(perlop) + numdatbns - (actfoinumsmp ./2)));
    nonacttimboiind = find(timboi <  (-minoffset + data.offset(perlop) + (actfoinumsmp ./ 2)) | timboi >= (-minoffset + data.offset(perlop) + numdatbns - (actfoinumsmp ./2)));
    acttimboi       = timboi(acttimboiind);
    numacttimboi    = length(acttimboi);
    if keep ==1
      cntpertoi(foilop,acttimboiind) = cntpertoi(foilop,acttimboiind) + 1;
    end
    for taplop = 1:numtap(foilop)
      if keep == 3
        cnt = taplop;
      elseif keep == 4
        % this once again assumes a fixed number of tapers per frequency
        cnt = (perlop-1)*numtap(1) + taplop;
      end
      autspctrmacttap = complex(zeros(numsgn,numacttimboi), zeros(numsgn,numacttimboi));
      if numacttimboi > 0
        for sgnlop = 1:numsgn
          dum = fftshift(ifft(datspctra(sgnlop,:) .* knlspctrmstr{foilop}(taplop,:),[],2));
          autspctrmacttap(sgnlop,:) = dum(acttimboi);
        end
      end
      if powflg
        powdum = 2.* abs(autspctrmacttap) .^ 2 ./ actfoinumsmp;
        if strcmp(cfg.taper, 'sine')
          powdum = powdum .* (1 - (((taplop - 1) ./ numtap(foilop)) .^ 2));
        end
        if keep == 1 && numacttimboi > 0
          powspctrm(:,foilop,acttimboiind) = powspctrm(:,foilop,acttimboiind) + reshape(powdum ./ numtap(foilop),[numsgn,1,numacttimboi]);
        elseif keep == 2 && numacttimboi > 0
          powspctrm(cnt,:,foilop,acttimboiind) = powspctrm(cnt,:,foilop,acttimboiind) + reshape(powdum ./ numtap(foilop),[1,numsgn,1,numacttimboi]);
          powspctrm(cnt,:,foilop,nonacttimboiind) = nan;
        elseif keep == 4 && numacttimboi > 0
          powspctrm(cnt,:,foilop,acttimboiind) = reshape(powdum,[1,numsgn,1,numacttimboi]);
          powspctrm(cnt,:,foilop,nonacttimboiind) = nan;
        elseif (keep == 4 || keep == 2) && numacttimboi == 0
          powspctrm(cnt,:,foilop,nonacttimboiind) = nan;
        end
      end
      if fftflg
        fourierdum = (autspctrmacttap) .* sqrt(2 ./ actfoinumsmp); %cf Numercial Receipes 13.4.9
        if keep == 1 && numacttimboi > 0
          fourierspctrm(:,foilop,acttimboiind) = fourierspctrm(:,foilop,acttimboiind) + reshape((fourierdum ./ numtap(foilop)),[numsgn,1,numacttimboi]);
        elseif keep == 2 && numacttimboi > 0
          fourierspctrm(cnt,:,foilop,acttimboiind) = fourierspctrm(cnt,:,foilop,acttimboiind) + reshape(fourierdum ./ numtap(foilop),[1,numsgn,1,numacttimboi]);
          fourierspctrm(cnt,:,foilop,nonacttimboiind) = nan;
        elseif keep == 4 && numacttimboi > 0
          fourierspctrm(cnt,:,foilop,acttimboiind) = reshape(fourierdum,[1,numsgn,1,numacttimboi]);
          fourierspctrm(cnt,:,foilop,nonacttimboiind) = nan;
        elseif (keep == 4 || keep == 2) && numacttimboi == 0
          fourierspctrm(cnt,:,foilop,nonacttimboiind) = nan;
        end
      end
      if csdflg
        csddum = 2.* (autspctrmacttap(cutdatindcmb(:,1),:) .* conj(autspctrmacttap(cutdatindcmb(:,2),:))) ./ actfoinumsmp;
        if keep == 1 && numacttimboi > 0
          crsspctrm(:,foilop,acttimboiind) = crsspctrm(:,foilop,acttimboiind) + reshape((csddum ./ numtap(foilop)),[numsgncmb,1,numacttimboi]);
        elseif keep == 2 && numacttimboi > 0
          crsspctrm(cnt,:,foilop,acttimboiind) = crsspctrm(cnt,:,foilop,acttimboiind) + reshape(csddum ./ numtap(foilop),[1,numsgncmb,1,numacttimboi]);
          crsspctrm(cnt,:,foilop,nonacttimboiind) = nan;
        elseif keep == 4 && numacttimboi > 0
          crsspctrm(cnt,:,foilop,acttimboiind) = reshape(csddum,[1,numsgncmb,1,numacttimboi]);
          crsspctrm(cnt,:,foilop,nonacttimboiind) = nan;
        elseif (keep == 4 || keep == 2) && numacttimboi == 0
          crsspctrm(cnt,:,foilop,nonacttimboiind) = nan;
        end
      end
    end % for taplop
    if calcdof
      dof(perlop,foilop,acttimboiind) = numtap(foilop);
    end
  end % for foilop
end % for perlop

if keep == 1
  warning off
  if powflg
    powspctrm(:,:,:) = powspctrm(:,:,:) ./ repmat(permute(cntpertoi,[3,1,2]),[numsgn,1,1]);
  end
  if fftflg
    fourierspctrm(:,:,:) = fourierspctrm(:,:,:) ./ repmat(permute(cntpertoi,[3,1,2]),[numsgn,1,1]);
  end
  if csdflg
    crsspctrm(:,:,:) = crsspctrm(:,:,:) ./ repmat(permute(cntpertoi,[3,1,2]),[numsgncmb,1,1]);
  end
  warning on
end

% collect the results
freq.label      = data.label(sgnindx);
freq.dimord     = dimord;
freq.freq       = cfg.foi;
freq.time       = toi;

if powflg
  freq.powspctrm  = powspctrm;
end

if csdflg
  freq.labelcmb   = cfg.channelcmb;
  freq.crsspctrm  = crsspctrm;
end

if fftflg
  freq.fourierspctrm  = fourierspctrm;
end

if calcdof
  freq.dof=2*dof;
end;

if keep == 2,
  freq.cumtapcnt = repmat(numtap(:)', [size(powspctrm,1) 1]);
elseif keep == 4,
  %all(numtap(1)==numtap)
  freq.cumtapcnt = repmat(numtap(1), [size(fourierspctrm,1)./numtap(1) 1]);
end

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
cfg.version.id = '$Id: freqanalysis_mtmconvol.m,v 1.44 2009/03/11 10:39:37 roboos Exp $';
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

