function [output] = freqdescriptives(cfg, freq)

% FREQDESCRIPTIVES computes descriptive statistics of the frequency
% or time-frequency decomposition of the EEG/MEG signal, such as
% the average power and the coherence.
%
% Use as
%   [freq] = freqdescriptives(cfg, freq)
%
% The data in freq should be organised in a structure as obtained from
% from the FREQANALYSIS function. The output structure is comparable
% to the input structure and can be used in most functions that require
% a freq input.
%
% The configuration options are
%   cfg.cohmethod     = 'coh' or 'plv' computes coherence or phase-locking-value (default = 'coh')
%   cfg.complex       = 'abs', 'complex', 'real', 'imag', 'absreal', 'absimag' or 'angle' (default = 'abs')
%   cfg.combinechan   = 'no' or 'planar' (default = 'no'), see below
%   cfg.combinemethod = 'svdfft', algorithm that is used to combine planar channels: the gradients
%                        are projected on the direction in which the average power over trials is maximal
%   cfg.variance      = 'yes' or 'no', estimate standard error in the standard way (default = 'no)
%                       works only for power
%   cfg.jackknife     = 'yes' or 'no', estimate standard error by means of the jack-knife (default = 'yes')
%                       for power and coherence
%   cfg.biascorrect   = 'yes' or 'no', calculate jackknife bias-corrected power and coherence (default = 'no')
%                       this option can only chosen if cfg.jackknife = 'yes'
%   cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'),
%                       see CHANNELSELECTION for details
%   cfg.channelcmb    = Mx2 cell-array with selection of channel pairs (default = {'all' 'all'}),
%                       see CHANNELCOMBINATION for details
%   cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.foilim        = [fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
%   cfg.toilim        = [tmin tmax] or 'all', to specify a subset of latencies (default = 'all')
%
% Coherence can only be computed if the input data contains either power-
% and cross-spectral densities, or the single-trial Fourier spectra. The
% option cfg.channelcmb only applies to Fourier data.
%
% A variance estimate can only be computed if results from trials and/or
% tapers have been kept.
%
% If the data consists of planar gradient data, you can specify
% cfg.combinechan = 'no'. In that case the power-spectra can be computed
% afterwards by using COMBINEPLANAR. However, if you want to compute
% coherence between two pairs of planar gradiometers, you should specify
% cfg.combinechan = 'planar'. This method only works if your input data
% contains Fourier spectra, i.e. you should specify cfg.output = 'fourier'
% in FREQANALYSIS.
%
% See also FREQANALYSIS, FREQSTATISTICS, FREQBASELINE

% FIXME: include pseudovalue-stuff
% FIXME: optional z-transformation of coherence: for this, mtmconvol has to output a cumtapcnt
% FIXME: foilim and toilim moeten op alle input data werken
%
% Undocumented local options:
% cfg.feedback
% cfg.keeptrials
% cfg.latency
% cfg.partchan = cell-array (default is empty)
% cfg.previous
% cfg.pseudovalue
% cfg.version

% Copyright (C) 2004-2006, Pascal Fries & Jan-Mathijs Schoffelen, F.C. Donders Centre
%
% $Log: freqdescriptives.m,v $
% Revision 1.59  2009/06/12 11:48:22  jansch
% added default for cfg.keepfourier
%
% Revision 1.58  2009/04/08 06:05:23  roboos
% give warning in case plv and input only one trial (solves the problem of Wendy)
%
% Revision 1.57  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.56  2009/01/12 13:05:20  sashae
% small change in call to checkconfig
%
% Revision 1.55  2008/11/28 17:33:19  sashae
% fixed bug in calculating sumcrsspctrm when input data has no time dim
%
% Revision 1.54  2008/11/27 08:48:39  kaigoe
% Speedup by removing for-loops. Added some comments regarding coherence.
% Further speedup possible, if SEM Jackknife estimation should not be
% computed (3x).
%
% Revision 1.52  2008/09/30 16:45:55  sashae
% checkconfig: checks if the input cfg is valid for this function
%
% Revision 1.51  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.50  2008/09/22 11:30:42  roboos
% added default keeptrials=no, which appears to have gone missing in the previous commit
%
% Revision 1.49  2008/09/10 14:11:56  jansch
% fixed an issue with keeptrials
%
% Revision 1.48  2008/05/15 08:18:55  roboos
% replaced strmatch with strcmp where applicable to fix bug in toilim='all' (thanks to Doug)
%
% Revision 1.47  2008/05/06 16:30:26  sashae
% change in trial selection, cfg.trials can be logical
%
% Revision 1.46  2008/01/31 09:42:37  roboos
% keep cumtapcnt if available
%
% Revision 1.45  2008/01/29 18:17:40  sashae
% added option for trial selection
% input data with dimord 'subj_' now also possible, treated similarly as 'rpt_'
%
% Revision 1.44  2007/09/05 09:46:47  jansch
% fixed bug in passing of feedback to fourier2crsspctrm
%
% Revision 1.43  2007/05/30 11:38:43  roboos
% implemented real and absreal for cfg.complex, can be used for instantaneous correlation
%
% Revision 1.42  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.41  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.40  2007/03/27 11:05:19  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.39  2007/02/27 13:34:55  roboos
% added some comments for ingnie to work on
%
% Revision 1.38  2006/11/29 10:07:56  erimar
% Removed crash-producing code that was a leftover of a copy-paste.
%
% Revision 1.37  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.36  2006/07/04 17:06:23  ingnie
% updated documentation
%
% Revision 1.35  2006/07/04 16:04:50  roboos
% renamed option 'jacknife' into 'jackknife' for consistency, maintain backward compatibility with cfgs and old data
%
% Revision 1.34  2006/06/22 07:09:00  erimar
% Corrected a bug involving missing brackets.
%
% Revision 1.33  2006/06/20 13:20:13  erimar
% Added the config-option 'biascorrect', which allows to calculate the
% jackknife bias-corrected power and coherence.
%
% Revision 1.32  2006/06/19 14:43:13  erimar
% Corrected some reshape-operations (by adding 1 as the last argument).
%
% Revision 1.31  2006/06/06 16:57:51  ingnie
% updated documentation
%
% Revision 1.30  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.29  2006/03/22 13:59:59  jansch
% implemented support to do svdfft on input-data containing fourier-spectra
%
% Revision 1.28  2006/03/20 11:23:10  jansch
% new implementation, merging PF's functionality with RO's functionality
%
% Revision 1.26  2006/02/27 11:22:45  roboos
% fixed incorrect variable name, data should be freq
%
% Revision 1.25  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.24  2006/02/09 10:06:27  roboos
% changed a help comment
%
% Revision 1.23  2006/02/01 12:26:00  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.22  2005/09/05 06:37:51  jansch
% small change in feedback
%
% Revision 1.21  2005/08/16 11:19:04  jansch
% moved some stuff in 'ro's implementation to a separate function
% fourier2crsspctrm to facilitate re-use

fieldtripdefs

% check if the input data is valid for this function
% freq = checkdata(freq, 'datatype', 'freq', 'feedback', 'yes', 'hascumtapcnt', 'no');
freq = checkdata(freq, 'datatype', 'freq', 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'renamed',     {'jacknife',   'jackknife'});

% determine some specific details of the input data
hascsd   = isfield(freq, 'crsspctrm');
hasdft   = isfield(freq, 'fourierspctrm');
haspow   = isfield(freq, 'powspctrm');
hasrpt   = ~isempty(strfind(freq.dimord, 'rpt')) || ~isempty(strfind(freq.dimord, 'subj'));
hastap   = ~isempty(strfind(freq.dimord, 'tap'));
hastim   = ~isempty(strfind(freq.dimord, 'time'));

% set the defaults
if ~isfield(cfg, 'feedback'),          cfg.feedback      = 'textbar';     end
if ~isfield(cfg, 'cohmethod'),         cfg.cohmethod     = 'coh';         end
if ~isfield(cfg, 'complex'),           cfg.complex       = 'abs';         end
if ~isfield(cfg, 'combinechan'),       cfg.combinechan   = 'no';          end
if ~isfield(cfg, 'combinemethod'),     cfg.combinemethod = 'svdfft';      end
if ~isfield(cfg, 'pseudovalue'),       cfg.pseudovalue   = 'no';          end
if ~isfield(cfg, 'variance'),          cfg.variance      = 'no';          end
if ~isfield(cfg, 'trials'),            cfg.trials        = 'all';         end
if ~isfield(cfg, 'channel'),           cfg.channel       = 'all';         end
if ~isfield(cfg, 'partchan'),          cfg.partchan      = {};            end
if ~isfield(cfg, 'foilim'),            cfg.foilim        = 'all';         end
if ~isfield(cfg, 'toilim'),            cfg.toilim        = 'all';         end
if ~isfield(cfg, 'keeptrials'),        cfg.keeptrials    = 'no';          end
if ~isfield(cfg, 'keepfourier'),       cfg.keepfourier   = 'no';          end

if ~isfield(cfg, 'channelcmb'),
  if hascsd
    cfg.channelcmb    = freq.labelcmb;
  elseif hasdft
    cfg.channelcmb    = {'all', 'all'};
  else
    cfg.channelcmb    = {};
  end
end
if ~isfield(cfg, 'jackknife')
  if hasrpt
    cfg.jackknife     = 'yes';
  else
    cfg.jackknife     = 'no';
  end
end
if ~isfield(cfg, 'biascorrect')
  cfg.biascorrect = 'no';
end

jckflg   = strcmp(cfg.jackknife,   'yes');
psdflg   = strcmp(cfg.pseudovalue, 'yes');
varflg   = strcmp(cfg.variance,    'yes');
plvflg   = strcmp(cfg.cohmethod,   'plv');
bcrflg   = strcmp(cfg.biascorrect, 'yes');

if sum([jckflg psdflg varflg]>1)
  error('you should specify only one of cfg.jackknife, cfg.pseudovalue or cfg.variance');
end

if ~hasrpt && (jckflg || psdflg || varflg),
  error('a variance-estimate or pseudovalue-estimate without repeated observations in the input is not possible');
end

if ~(hascsd || hasdft) && ~isempty(cfg.partchan),
  error('cannot do partialisation without cross-spectra or fourier-spectra in the input data');
end

if strcmp(cfg.combinechan, 'planar') && ~hasdft,
  error('cfg.combinechan=''planar'' requires Fourier spectra in the input data');
end

if ~jckflg && bcrflg
  error('You can only calculate jackknife bias-corrected power and coherence estimates if cfg.jackknife=''yes''.');
end

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  if hastap
    error('trial selection is not possible for input data with a ''tap'' dimension');
  elseif ~hasrpt
    error('trial selection requires input data with repeated observations');
  end
  if hascsd
    if islogical(cfg.trials),  cfg.trials=find(cfg.trials);  end
    fprintf('selecting %d trials\n', length(cfg.trials));
    freq.crsspctrm=freq.crsspctrm(cfg.trials,:,:,:);
  end
  if haspow
    if islogical(cfg.trials),  cfg.trials=find(cfg.trials);  end
    fprintf('selecting %d trials\n', length(cfg.trials));
    freq.powspctrm=freq.powspctrm(cfg.trials,:,:,:);
  end
  % update the trial definition (trl)
  if isfield(freq, 'cfg') % try to locate the trl in the nested configuration
    trl = findcfg(freq.cfg, 'trl');
  else
    trl = [];
  end
  if isempty(trl)
    % a trial definition is expected in each continuous data set
    warning('could not locate the trial definition ''trl'' in the data structure');
  else
    cfg.trlold=trl;
    cfg.trl=trl(cfg.trials,:);
  end
end

%FIXME: use tokenize
if hasrpt,
  if strfind(freq.dimord, 'rpttap'),
    newdimord = freq.dimord(8:end);
  elseif strfind(freq.dimord, 'rpt'),
    newdimord = freq.dimord(5:end);
  elseif strfind(freq.dimord, 'subj'),
    newdimord = freq.dimord(6:end);
  end
else
  if hascsd
    freq.crsspctrm = reshape(freq.crsspctrm, [1 size(freq.crsspctrm)]);
  end
  if haspow
    freq.powspctrm = reshape(freq.powspctrm, [1 size(freq.powspctrm)]);
  end
  newdimord = freq.dimord;  % without the rpt
end

%check which frequency bins are requested
if strcmp(cfg.foilim, 'all'),
  fbin = [1:length(freq.freq)];
else
  fbin = [nearest(freq.freq, cfg.foilim(1)):nearest(freq.freq, cfg.foilim(2))];
end

%check which time bins are requested
if ~hastim,
  tbin   = 1;
elseif hastim && strcmp(cfg.toilim, 'all'),
  tbin   = [1:length(freq.time)];
else
  tbin   = [nearest(freq.time, cfg.latency(1)):nearest(freq.time, cfg.latency(2))];
end
hastim  = length(tbin)>1;
Nfrq    = length(fbin);
Ntim    = length(tbin);

% preprocessing of Fourier spectra
if hasdft,
  Nrpttap = size(freq.fourierspctrm,1);
  if ~strcmp(cfg.combinechan, 'no')
    if ~strcmp(cfg.combinemethod, 'svdfft'),
      error(sprintf('cfg.combinemethod=''%s'' is not yet implemented', cfg.combinemethod));
    end
    % inputlabel2outputlabel will check whether combinechan = 'planar' and act accordingly
    tmpcfg = [];
    tmpcfg.combinechan = cfg.combinechan;
    [outputlabel, outputindex] = inputlabel2outputlabel(tmpcfg, freq);
    fourier = zeros(Nrpttap, length(outputindex), Nfrq, Ntim);
    progress('init', cfg.feedback, 'computing svdfft for combined channels');
    for j = 1:length(outputindex)
      progress(j/length(outputindex));
      if length(outputindex{j}) > 1,
        for k = 1:Nfrq
          for m = 1:Ntim
            tmpsel = find(~isnan(freq.fourierspctrm(:, outputindex{j}(1), fbin(k), tbin(m))));
            fourier(tmpsel, j, k, m) = svdfft(transpose(freq.fourierspctrm(tmpsel, outputindex{j}, fbin(k), tbin(m))), 1);
            fourier(setdiff(1:size(freq.fourierspctrm,1),tmpsel), j, k, m) = nan;
          end
        end
      else
        fourier(:, j, :, :) = freq.fourierspctrm(:, outputindex{j}, fbin, tbin);
      end
    end
    progress('close');
    freq.fourierspctrm = fourier;
    freq.label         = outputlabel;
    try, freq.time     = freq.time(tbin); end
    freq.freq          = freq.freq(fbin);
    clear fourier;
  else
    % FIXME do tbin and fbin selection, also when combinechan=no
  end
end

if haspow,
  % FIXME, do fbin and tbin selection
end

if hascsd,
  % FIXME, do fbin and tbin selection
end

%channels should be specified as channels in the desired output
cfg.channel    = channelselection(cfg.channel,      freq.label);
cfg.partchan   = channelselection(cfg.partchan,     freq.label);

cfg.channelcmb = channelcombination(cfg.channelcmb, freq.label);
%further preprocessing of fourierspectra, not needed to consider planar combinations anymore
if hasdft,
  %convert to cross-spectra
  partindx = match_str(freq.label, cfg.partchan);
  chnindx  = match_str(freq.label, cfg.channel);
  %check which channels are needed
  chancmb  = cfg.channelcmb;
  for j = 1:length(cfg.partchan),
    chancmb       = [cfg.channel repmat(cfg.partchan(j), [length(cfg.channel) 1]); chancmb];
  end
  chnindx = unique([chnindx; partindx]);
  tmpcfg = [];
  tmpcfg.channel    = freq.label(chnindx);
  tmpcfg.channelcmb = chancmb;
  tmpcfg.keeptrials = 'yes';
  tmpcfg.feedback   = cfg.feedback;
  freq   = fourier2crsspctrm(tmpcfg, freq);
  hascsd = isfield(freq, 'crsspctrm');
  haspow = isfield(freq, 'powspctrm');
end

%concatenate cross-spectra and power-spectra
if hascsd && haspow,
  freq.crsspctrm = cat(2,freq.powspctrm,freq.crsspctrm);
  freq.labelcmb  = [freq.label(:) freq.label(:); freq.labelcmb];
  freq = rmfield(freq,'powspctrm');
elseif haspow,
  freq.crsspctrm = freq.powspctrm;
  freq.labelcmb  = [freq.label(:) freq.label(:)];
  freq = rmfield(freq,'powspctrm');
end

%normalise cross-spectra
if plvflg,
  if size(freq.crsspctrm,1)==1
    warning('It seems that the data only contains a single trial or an average, which makes it impossible to compute the phase-locking value. Use cfg.keeptrials=''yes'' in freqanalysis.')
  end
  freq.crsspctrm = freq.crsspctrm./abs(freq.crsspctrm);
end

%prepare some stuff for jack-knifing
Nrpt = size(freq.crsspctrm,1);
Ncmb = size(freq.crsspctrm,2);
ind  = isfinite(freq.crsspctrm);
dof  = sum(ind,1);
if jckflg || psdflg || varflg
  % these can only be computed if there are enough trials
  dof(find(dof<3)) = nan;
end

% new efficient version
sumcrsspctrm(1,1:Ncmb,1:Nfrq,1:Ntim) = nansum(freq.crsspctrm,1);
%% does the same as the old inefficient version
% sumcrsspctrm = complex(zeros(1,Ncmb,Nfrq,Ntim), ...
%                         zeros(1,Ncmb,Nfrq,Ntim));
% for j = 1:Ncmb
%   sumcrsspctrm(1,j,:,:) = nansum(freq.crsspctrm(:,j,:,:),1);
% end
avgcrsspctrm = sumcrsspctrm./dof;

%compute leave-one-out averages
if jckflg || psdflg,
  progress('init', cfg.feedback, 'computing the leave-one-out averages');
  for j = 1:Nrpt
    progress(j/Nrpt);
    freq.crsspctrm(j,:,:,:) = (sumcrsspctrm - freq.crsspctrm(j,:,:,:))./(dof-double(ind(j,:,:,:)));
  end
  progress('close');
elseif varflg,
  %do nothing
  warning('variance computation only interpretable for power estimate');
end

%data-dimensionality
Nrpt = size(freq.crsspctrm,1);
Nsgn = length(freq.label);
Nfrq = size(freq.crsspctrm,3);
Ntim = size(freq.crsspctrm,4);

%partialise the cross-spectra for the leave-one out averages
if ~isempty(cfg.partchan) && (jckflg || psdflg),
  [freq] = partialisation(cfg, freq, 'crsspctrm');
elseif ~isempty(cfg.partchan) && (varflg)
  error('you cannot partialise out channels without leave-one-out averaging');
end

%partialise the cross-spectra for the average
if ~isempty(cfg.partchan),
  tmpfreq           = freq;
  tmpfreq.crsspctrm = avgcrsspctrm;
  [tmpfreq]         = partialisation(cfg, tmpfreq, 'crsspctrm');
  avgcrsspctrm      = tmpfreq.crsspctrm;
  clear tmpfreq;
end

%split cross-spectra into true cross-spectra and autospectra

% determine the corresponding indices of all channel combinations
cmbindx = zeros(size(freq.labelcmb));
for k=1:size(freq.labelcmb,1)
  cmbindx(k,1) = find(strcmp(freq.labelcmb(k,1), freq.label));
  cmbindx(k,2) = find(strcmp(freq.labelcmb(k,2), freq.label));
  % this works the same, but is much slower in Octave
  % sgncmbindx(k,1) = strmatch(freq.labelcmb(k,1), freq.label, 'exact');
  % sgncmbindx(k,2) = strmatch(freq.labelcmb(k,2), freq.label, 'exact');
end
autocmb = cmbindx(:,1) == cmbindx(:,2);
autocmb = find(autocmb);

avgpowspctrm   = avgcrsspctrm(:, autocmb, :, :);
avgcrsspctrm   = avgcrsspctrm(:, setdiff(1:size(cmbindx,1), autocmb), :, :);
freq.powspctrm = freq.crsspctrm(:, autocmb, :, :);
freq.crsspctrm = freq.crsspctrm(:, setdiff(1:size(cmbindx,1), autocmb), :, :);
freq.label     = freq.label(cmbindx(autocmb,1));
freq.labelcmb  = freq.labelcmb(setdiff(1:size(cmbindx,1), autocmb), :);
dofcsd         = dof(:, setdiff(1:size(cmbindx,1), autocmb), :, :);
dofpow         = dof(:, autocmb, :, :);
haspow         = 1;  % from now on it again contains a seperate power spectrum

%compute the coherence
if hascsd,
  % determine the corresponding indices of all channel combinations
  cmbindx = zeros(size(cfg.channelcmb));
  for k=1:size(cfg.channelcmb,1)
    cmbindx(k,1) = find(strcmp(cfg.channelcmb(k,1), freq.label));
    cmbindx(k,2) = find(strcmp(cfg.channelcmb(k,2), freq.label));
    % this works the same, but is much slower in Octave
    % sgncmbindx(k,1) = strmatch(cfg.channelcmb(k,1), data.label, 'exact');
    % sgncmbindx(k,2) = strmatch(cfg.channelcmb(k,2), data.label, 'exact');
  end

  % TODO: THIS is the COHERENCE
  % AT LEAST THIS IS WHAT LINE
  %  621: output.cohspctrm = reshape(avgcrsspctrm, [sizcrs(2:end), 1]);
  % says

  fprintf('computing the coherence\n');
  % new efficient version
  avgcrsspctrm = avgcrsspctrm./sqrt( avgpowspctrm(:, cmbindx(:,1), :, :) .* avgpowspctrm(:, cmbindx(:,2), :, :) );
  % does the same as the old, unefficient version:
  %   for j = 1 %average consists of 1 repetition
  %     for k = 1:Nfrq
  %       for m = 1:Ntim
  %         avgcrsspctrm(j, :, k, m) = avgcrsspctrm(j, :, k, m)./sqrt( avgpowspctrm(j, cmbindx(:,1), k, m) .* avgpowspctrm(j, cmbindx(:,2), k, m) );
  %       end
  %     end
  %   end

  % TODO: If I see this correctly, thats not the coherence. Instead, the
  % avgcrsspctrm computed above seems to be the coherence...
  % TODO: It this is true, than avoid computing it, if we do not need it
  % (takes quite long)

  fprintf('computing the coherence for SEM\n');
  % new efficient version
  freq.crsspctrm = freq.crsspctrm ./sqrt( freq.powspctrm(:, cmbindx(:,1), :, :) .* freq.powspctrm(:, cmbindx(:,2), :, :) );
  % performs the same operation as the old but inefficient version
  % (for larger datassets, >50x faster - Order: O(>n))
  %   for j = 1:Nrpt
  %     progress(j/Nrpt);
  %     for k = 1:Nfrq
  %       for m = 1:Ntim
  %         freq.crsspctrm(j, :, k, m) = freq.crsspctrm(j, :, k, m)./sqrt( freq.powspctrm(j, cmbindx(:,1), k, m) .* freq.powspctrm(j, cmbindx(:,2), k, m) );
  %       end
  %     end
  %   end


  %put output coherence in correct format
  switch cfg.complex
    case 'complex'
      % leave all values complex
    case 'abs'
      % convert to absolute values
      freq.crsspctrm = abs(freq.crsspctrm);
      avgcrsspctrm   = abs(avgcrsspctrm);
    case 'real'
      freq.crsspctrm = real(freq.crsspctrm);
      avgcrsspctrm   = real(avgcrsspctrm);
    case 'imag'
      freq.crsspctrm = imag(freq.crsspctrm);
      avgcrsspctrm   = imag(avgcrsspctrm);
    case 'absreal'
      freq.crsspctrm = abs(real(freq.crsspctrm));
      avgcrsspctrm   = abs(real(avgcrsspctrm));
    case 'absimag'
      freq.crsspctrm = abs(imag(freq.crsspctrm));
      avgcrsspctrm   = abs(imag(avgcrsspctrm));
    case 'angle'
      freq.crsspctrm = angle(freq.crsspctrm);
      avgcrsspctrm   = angle(avgcrsspctrm);
    otherwise
      error(sprintf('method ''%s'' is not implemented for cfg.complex', cfg.complex));
  end
end

%compute sem
if haspow,
  powspctrmsem = zeros(size(freq.powspctrm,2),Nfrq,Ntim);
  if bcrflg
    powspctrmbcr = zeros(1,size(freq.powspctrm,2),Nfrq,Ntim);
  end
  if jckflg,
    for j = 1:Ntim
      powspctrmsem(:,:,j) = squeeze(nanstd(freq.powspctrm(:, :, :, j), 1, 1).*sqrt(dofpow(:,:,:,j)-1)); %cf Efron p.141
      if bcrflg
        % temporarily store the average over the jackknife samples in powspctrmbcr
        powspctrmbcr(:,:,:,j) = nanmean(freq.powspctrm(:, :, :, j), 1);
      end;
    end
  elseif varflg,
    for j = 1:Ntim
      powspctrmsem(:,:,j) = squeeze(nanstd(freq.powspctrm(:, :, :, j), 0, 1)./sqrt(dofpow(:,:,:,j)));
    end
  end
end

if hascsd,
  if jckflg
    cohspctrmsem = zeros(size(freq.crsspctrm,2),Nfrq,Ntim);
    if bcrflg
      cohspctrmbcr = zeros(1,size(freq.crsspctrm,2),Nfrq,Ntim);
    end
    for j = 1:Ntim
      cohspctrmsem(:,:,j) = squeeze(nanstd(freq.crsspctrm(:, :, :, j), 1, 1).*sqrt(dofcsd(:,:,:,j)-1));
      if bcrflg
        % temporarily store the average over the jackknife samples in cohspctrmbcr
        cohspctrmbcr(:,:,:,j) = nanmean(freq.crsspctrm(:, :, :, j), 1);
      end;
    end
  elseif varflg,
    fprintf('cannot compute the variance of the coherence without leave-one-out resampling\n');
  end
end

%create the output-structure
output                 = [];
output.dimord          = newdimord;
output.freq            = freq.freq;
if hastim,
  output.time = freq.time;
end
output.label           = freq.label;

if strcmp(cfg.keeptrials, 'no'),
  sizpow                 = size(avgpowspctrm);
  output.powspctrm       = reshape(avgpowspctrm, [sizpow(2:end), 1]);
  if jckflg || varflg,
    output.powspctrmsem = powspctrmsem;
  end
  if hascsd,
    output.labelcmb     = freq.labelcmb;
    sizcrs = size(avgcrsspctrm);
    if plvflg,
      % rename to plv
      output.plvspctrm = reshape(avgcrsspctrm, [sizcrs(2:end), 1]);
    else
      output.cohspctrm = reshape(avgcrsspctrm, [sizcrs(2:end), 1]);
    end
    if jckflg && plvflg,
      % rename to plv
      output.plvspctrmsem = cohspctrmsem;
    elseif jckflg
      output.cohspctrmsem = cohspctrmsem;
    end
  end
  if hascsd,
    output.dof = reshape(dofcsd, [sizcrs(2:end), 1]);
  else
    output.dof = reshape(dofpow, [sizpow(2:end), 1]);
  end
else
  output.powspctrm = freq.powspctrm;
  output.dimord    = ['rpt_',output.dimord];
  if hascsd,
    output.crsspctrm = freq.crsspctrm;
    output.labelcmb  = freq.labelcmb;
    output.cumtapcnt = freq.cumtapcnt;
  end
  if strcmp(cfg.keepfourier, 'yes'),
    output.fourierspctrm = freq.fourierspctrm;
    output.dimord        = ['rpttap_',output.dimord(5:end)];
    output.cumtapcnt     = freq.cumtapcnt;
    output               = rmfield(output, 'powspctrm');
  end
end
try, output.grad       = freq.grad; end

%FIXME: pseudovalue
%if psdflg,
%  Nrpt = size(pow,1) - 1;
%  output.pseudo.powspctrm = repmat(Nrpt.*pow(end,:,:,:), [Nrpt 1 1 1]) - (Nrpt-1).*pow(:,:,:,:);
%end
%
%
%  if psdflg,
%    output.pseudo.cohspctrm = repmat(Nrpt.*coh, [Nrpt 1 1 1]) - (Nrpt-1).*jckcoh;
%  end

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
cfg.version.id = '$Id: freqdescriptives.m,v 1.59 2009/06/12 11:48:22 jansch Exp $';
try, cfg.previous = freq.cfg; end

% remember the configuration details
output.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that partialises cross-spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [freq] = partialisation(cfg, freq, field)

crsspctrm = getfield(freq, field);
Nrpt      = size(crsspctrm,1);
Nsgn      = size(crsspctrm,2);
Nfrq      = size(crsspctrm,3);
Ntim      = size(crsspctrm,4);
%create signalcombination-matrix, each non-nan entry corresponding
%with the index of the signalcombination in the cross-spectrum
%negative entries should be conjugated before further processing.
%the underlying idea is to temporarily convert the cross-spectra into
%a matrix-shape, so that the partialisation can take place efficiently
%this mainly involves a lot of careful book-keeping
fprintf('indexing the channels and channel combinations\n');
sgnindx    = [1:length(freq.label)]';
cmbindx = zeros(size(freq.labelcmb));
for i=1:size(freq.labelcmb,1)
  cmbindx(i,1) = find(strcmp(freq.labelcmb(i,1), freq.label));
  cmbindx(i,2) = find(strcmp(freq.labelcmb(i,2), freq.label));
  % this works the same, but is much slower in Octave
  % cmbindx(i,1) = strmatch(freq.labelcmb{i,1}, freq.label, 'exact');
  % cmbindx(i,2) = strmatch(freq.labelcmb{i,2}, freq.label, 'exact');
end
sgncmbmat  = nan*zeros(length(sgnindx),length(sgnindx));
for i=1:size(cmbindx,1)
  sgncmbmat(cmbindx(i,2),cmbindx(i,1)) = -i;
  sgncmbmat(cmbindx(i,1),cmbindx(i,2)) = i;
end

fprintf('checking whether partialisation of the requested channel(s) is possible\n');
prtindx = [];
for j = 1:length(cfg.partchan)
  prtindx = [prtindx; strmatch(cfg.partchan{j}, freq.label, 'exact')];
end
if isempty(prtindx), warning('no partialisation will be performed, requested channel(s) not present in data'); end;

prtindx     = sort(prtindx);
sgnindx     = setdiff(sgnindx, prtindx);
prtsgnmat   = sgncmbmat(prtindx,prtindx);
rstsgnmat   = sgncmbmat(sgnindx,sgnindx);
Nprtsgn     = length(prtindx);
Nrstsgn     = length(sgnindx);
prtsgnlabel = freq.label(prtindx);
label       = freq.label(sgnindx);

%the reformatting of the cross-spectral densities leads to a re-ordering of the signal-combinations
[ind1, ind2] = find(~isnan(rstsgnmat));
labelcmb     = [label(ind1) label(ind2)];
cmbindx      = [ind1 ind2];
for j = 1:size(cmbindx,1)
  duplicate(j,1) = find([cmbindx(:,1)==cmbindx(j,2)] .* [cmbindx(:,2)==cmbindx(j,1)]);
end %but we end up with a redundant amount of cross-spectra
powindx = find(duplicate==duplicate(duplicate)); %indices of auto-spectra
sel     = setdiff(1:size(cmbindx,1),powindx);
duplicate(powindx) = nan;
for j = 1:length(duplicate)
  if ~isnan(duplicate(j)), duplicate(find(duplicate==j)) = nan; end
end

%these cross-spectral density matrices should be complete
if any(isnan(prtsgnmat)),    error('partialisation of the requested channel(s) is not possible'); end;
prtrstsgnmat = sgncmbmat(prtindx,sgnindx);
if any(isnan(prtrstsgnmat)), error('partialisation of the requested channel(s) is not possible'); end;
rstprtsgnmat = sgncmbmat(sgnindx,prtindx);
if any(isnan(rstprtsgnmat)), error('partialisation of the requested channel(s) is not possible'); end;

%allocate memory
csdrr = nan+complex(zeros(Nrstsgn.^2,Nfrq,Ntim)      , zeros(Nrstsgn.^2,Nfrq,Ntim));
csdpp = nan+complex(zeros(Nprtsgn.^2,Nfrq,Ntim)      , zeros(Nprtsgn.^2,Nfrq,Ntim));
csdrp = nan+complex(zeros(Nrstsgn.*Nprtsgn,Nfrq,Ntim), zeros(Nrstsgn.*Nprtsgn,Nfrq,Ntim));
csdpr = nan+complex(zeros(Nprtsgn.*Nrstsgn,Nfrq,Ntim), zeros(Nprtsgn.*Nrstsgn,Nfrq,Ntim));

progress('init', cfg.feedback, 'partialising out the requested channel(s)')
for m = 1:size(crsspctrm,1)
  progress(m/size(crsspctrm,1));
  %a bit of hocus-pocus with the cross-spectral densities
  csdrr(find(~isnan(rstsgnmat)),:,:)  = squeeze(crsspctrm(m,abs(rstsgnmat(find(~isnan(rstsgnmat)))),:,:));
  csdpp(:) = squeeze(crsspctrm(m,abs(prtsgnmat),:,:));
  csdrp(:) = squeeze(crsspctrm(m,abs(rstprtsgnmat),:,:));
  csdpr(:) = squeeze(crsspctrm(m,abs(prtrstsgnmat),:,:));

  %take care of the conjugates
  csdrr(find(rstsgnmat<0),:,:)     = conj(csdrr(find(rstsgnmat<0),:,:));
  csdpp(find(prtsgnmat<0),:,:)     = conj(csdpp(find(prtsgnmat<0),:,:));
  csdrp(find(rstprtsgnmat<0),:,:)  = conj(csdrp(find(rstprtsgnmat<0),:,:));
  csdpr(find(prtrstsgnmat<0),:,:)  = conj(csdpr(find(prtrstsgnmat<0),:,:));

  %perform the partialisation
  for j = 1:Nfrq
    for k = 1:Ntim
      rr  = reshape(squeeze(csdrr(:,j,k)),[Nrstsgn Nrstsgn]);
      rp  = reshape(squeeze(csdrp(:,j,k)),[Nrstsgn Nprtsgn]);
      pp  = reshape(squeeze(csdpp(:,j,k)),[Nprtsgn Nprtsgn]);
      pr  = reshape(squeeze(csdpr(:,j,k)),[Nprtsgn Nrstsgn]);
      rrp = rr - rp*pinv(pp)*pr;
      csdp(m,:,j,k) = rrp(find(~isnan(rrp)));
    end
  end
end
progress('close');

freq.label     = label;
freq.labelcmb  = labelcmb([powindx; find(~isnan(duplicate))], :);
freq           = setfield(freq, field, csdp(:, [powindx; find(~isnan(duplicate))], :, :));
