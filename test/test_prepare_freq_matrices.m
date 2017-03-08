function test_prepare_freq_matrices

% WALLTIME 00:10:00
% MEM 1000mb

% TEST prepare_freq_matrices ft_sourceanalysis

datadir = dccnpath('/home/common/matlab/fieldtrip/data/test/latest/freq/meg');

cd(dccnpath('/home/common/matlab/fieldtrip/private'));

% fourier data, multiple trials
load(fullfile(datadir,'freq_mtmfft_fourier_trl_ctf275.mat'));

cfg           = [];
cfg.frequency = 5;
cfg.channel   = ft_channelselection('MEG',freq.label);
[a1,a2,a3,a4] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));

cfg.frequency = 5.4;
[a1,a2,a3,a4,cfg1] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4,cfg2] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));

cfg.frequency = 10;
cfg.refchan   = 'BR1';
[a1,a2,a3,a4,cfg1] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4,cfg2] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));
assert(isequal(a2,b2) ||  norm(a2-b2)<eps^2);
assert(isequal(a3,b3) ||  norm(a3-b3)<eps^2);

cfg.frequency = 10.6;
cfg.refchan   = 'BR1';
[a1,a2,a3,a4,cfg1] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4,cfg2] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));
assert(isequal(a2,b2) ||  norm(a2-b2)<eps^2);
assert(isequal(a3,b3) ||  norm(a3-b3)<eps^2);

% powandcsd data, multiple trials
load(fullfile(datadir,'freq_mtmfft_powandcsd_trl_ctf275.mat'));

cfg           = [];
cfg.frequency = 5;
cfg.channel   = ft_channelselection('MEG',freq.label);
[a1,a2,a3,a4] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));

cfg.frequency = 5.4;
[a1,a2,a3,a4,cfg1] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4,cfg2] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));

cfg.frequency = 10;
cfg.refchan   = 'BR1';
[a1,a2,a3,a4,cfg1] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4,cfg2] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));
assert(isequal(a2,b2));
assert(isequal(a3,b3));

cfg.frequency = 10.6;
cfg.refchan   = 'BR1';
[a1,a2,a3,a4,cfg1] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4,cfg2] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));
assert(isequal(a2,b2));
assert(isequal(a3,b3));

% powandcsd data, multiple trials and time
load(fullfile(datadir,'freq_mtmconvol_powandcsd_trl_ctf275.mat'));

cfg           = [];
cfg.frequency = 6;
cfg.latency   = 0.5;
cfg.channel   = ft_channelselection('MEG',freq.label);
[a1,a2,a3,a4] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));

cfg.frequency = 5.5;
[a1,a2,a3,a4,cfg1] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4,cfg2] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));

cfg.frequency = 6;
cfg.latency   = 0.54;
[a1,a2,a3,a4,cfg1] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4,cfg2] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));


cfg.frequency = 10;
cfg.refchan   = 'BR1';
[a1,a2,a3,a4,cfg1] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4,cfg2] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));
assert(isequal(a2,b2));
assert(isequal(a3,b3));

cfg.frequency = 10.6;
cfg.refchan   = 'BR1';
[a1,a2,a3,a4,cfg1] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4,cfg2] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));
assert(isequal(a2,b2));
assert(isequal(a3,b3));

% fourier data, multiple trials and time
load(fullfile(datadir,'freq_mtmconvol_fourier_trl_ctf275.mat'));

cfg           = [];
cfg.frequency = 6;
cfg.latency   = 0.5;
cfg.channel   = ft_channelselection('MEG',freq.label);
[a1,a2,a3,a4] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));

cfg.frequency = 5.5;
[a1,a2,a3,a4,cfg1] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4,cfg2] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));

cfg.frequency = 6;
cfg.latency   = 0.54;
[a1,a2,a3,a4,cfg1] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4,cfg2] = prepare_freq_matrices_old(cfg, freq);
assert(isequal(a1,b1));


cfg.frequency = 10;
cfg.refchan   = 'BR1';
[a1,a2,a3,a4,cfg1] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4,cfg2] = prepare_freq_matrices_old(cfg, freq);
assert(isalmostequal(a1,b1,'reltol',1e-9));
assert(isalmostequal(a2,b2,'reltol',1e-9));
assert(isalmostequal(a3,b3,'reltol',1e-9));

cfg.frequency = 10.6;
cfg.refchan   = 'BR1';
[a1,a2,a3,a4,cfg1] = prepare_freq_matrices(cfg, freq);
[b1,b2,b3,b4,cfg2] = prepare_freq_matrices_old(cfg, freq);
assert(isalmostequal(a1,b1,'reltol',1e-9));
assert(isalmostequal(a2,b2,'reltol',1e-9));
assert(isalmostequal(a3,b3,'reltol',1e-9));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BELOW IS THE OLD CODE

function [Cf, Cr, Pr, Ntrials, cfg] = prepare_freq_matrices_old(cfg, freq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that converts a freq structure into Cf, Cr and Pr
% this is used in sourecanalysis
%
% This function returns data matrices with a channel order that is consistent
% with the original channel order in the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the defaults
if ~isfield(cfg, 'dicsfix'), cfg.dicsfix = 'yes'; end
if ~isfield(cfg, 'quickflag'), cfg.quickflag = 0; end
if ~isfield(cfg, 'refchan'), cfg.refchan = []; end

quickflag = cfg.quickflag==1;

Cf = [];
Cr = [];
Pr = [];

% select the latency of interest for time-frequency data
if strcmp(freq.dimord, 'chan_freq_time')
  tbin = nearest(freq.time, cfg.latency);
  fprintf('selecting timeslice %d\n', tbin);
  freq.time = freq.time(tbin);
  % remove all other latencies from the data structure and reduce the number of dimensions
  if isfield(freq, 'powspctrm'),     freq.powspctrm     = squeeze(freq.powspctrm(:,:,tbin));     end;
  if isfield(freq, 'crsspctrm'),     freq.crsspctrm     = squeeze(freq.crsspctrm(:,:,tbin));     end;
  if isfield(freq, 'fourierspctrm'), freq.fourierspctrm = squeeze(freq.fourierspctrm(:,:,tbin)); end;
  freq.dimord = freq.dimord(1:(end-5));  % remove the '_time' part
elseif strcmp(freq.dimord, 'rpt_chan_freq_time') || strcmp(freq.dimord, 'rpttap_chan_freq_time')
  tbin = nearest(freq.time, cfg.latency);
  fprintf('selecting timeslice %d\n', tbin);
  freq.time = freq.time(tbin);
  % remove all other latencies from the data structure and reduce the number of dimensions
  if isfield(freq, 'powspctrm'),    freq.powspctrm     = squeeze(freq.powspctrm(:,:,:,tbin));      end;
  if isfield(freq, 'crsspctrm'),    freq.crsspctrm     = squeeze(freq.crsspctrm(:,:,:,tbin));      end;
  if isfield(freq, 'fourierspctrm') freq.fourierspctrm = squeeze(freq.fourierspctrm(:,:,:,tbin));  end;
  freq.dimord = freq.dimord(1:(end-5));  % remove the '_time' part
else
  tbin = [];
end

% the time-frequency latency has already been squeezed away (see above)
if strcmp(freq.dimord, 'chan_freq') || strcmp(freq.dimord, 'chancmb_freq') || strcmp(freq.dimord, 'chan_chan_freq') || strcmp(freq.dimord, 'chan_chan_freq_time')
  Ntrials = 1;
elseif strcmp(freq.dimord, 'rpt_chan_freq') || strcmp(freq.dimord, 'rpt_chancmb_freq') || strcmp(freq.dimord, 'rpt_chan_chan_freq')
  Ntrials = size(freq.cumtapcnt,1);
elseif strcmp(freq.dimord, 'rpttap_chan_freq') || strcmp(freq.dimord, 'rpttap_chancmb_freq') || strcmp(freq.dimord, 'rpttap_chan_chan_freq')
  Ntrials = size(freq.cumtapcnt,1);
elseif strcmp(freq.dimord, 'rpttap_chan_freq_time') || strcmp(freq.dimord, 'rpttap_chancmb_freq_time') || strcmp(freq.dimord, 'rpttap_chan_chan_freq_time')
  Ntrials = size(freq.cumtapcnt,1);
else
  error('unrecognized dimord for frequency data');
end

% find the frequency of interest
fbin = nearest(freq.freq, cfg.frequency);

if isfield(freq, 'powspctrm') && isfield(freq, 'crsspctrm')
  % use the power and cross spectrum and construct a square matrix
  
  % find the index of each sensor channel into powspctrm
  % keep the channel order of the cfg
  [dum, powspctrmindx] = match_str(cfg.channel, freq.label);
  Nchans = length(powspctrmindx);
  
  % find the index of each sensor channel combination into crsspctrm
  % keep the channel order of the cfg
  crsspctrmindx = zeros(Nchans);
  for sgncmblop=1:size(freq.labelcmb,1)
    ch1 = find(strcmp(cfg.channel, freq.labelcmb(sgncmblop,1)));
    ch2 = find(strcmp(cfg.channel, freq.labelcmb(sgncmblop,2)));
    if ~isempty(ch1) && ~isempty(ch2)
      % this square matrix contains the indices into the signal combinations
      crsspctrmindx(ch1,ch2) = sgncmblop;
    end
  end
  
  % this complex rearrangement of channel indices transforms the CSDs into a square matrix
  if strcmp(freq.dimord, 'chan_freq') || strcmp(freq.dimord, 'chancmb_freq')
    % FIXME this fails in case dimord=rpt_chan_freq and only 1 trial
    Cf = complex(nan(Nchans,Nchans));
    % first use the complex conjugate for all reversed signal combinations
    Cf(find(crsspctrmindx)) = freq.crsspctrm(crsspctrmindx(find(crsspctrmindx)), fbin);
    Cf = ctranspose(Cf);
    % and then get get the csd for all signal combinations
    Cf(find(crsspctrmindx)) = freq.crsspctrm(crsspctrmindx(find(crsspctrmindx)), fbin);
    % put the power on the diagonal
    Cf(find(eye(Nchans))) = freq.powspctrm(powspctrmindx, fbin);
  else
    Cf  = complex(nan(Ntrials,Nchans,Nchans));
    tmp = complex(nan(Nchans,Nchans));
    for trial=1:Ntrials
      % first use the complex conjugate for all signal combinations reversed
      tmp(find(crsspctrmindx)) = freq.crsspctrm(trial, crsspctrmindx(find(crsspctrmindx)), fbin);
      tmp = ctranspose(tmp);
      % and then get get the csd for all signal combinations
      tmp(find(crsspctrmindx)) = freq.crsspctrm(trial, crsspctrmindx(find(crsspctrmindx)), fbin);
      % put the power on the diagonal
      tmp(find(eye(Nchans))) = freq.powspctrm(trial, powspctrmindx, fbin);
      Cf(trial,:,:) = tmp;
    end
  end
  
  % do a sanity check on the cross-spectral-density matrix
  if any(isnan(Cf(:)))
    error('The cross-spectral-density matrix is not complete');
  end
  
  if isfield(cfg, 'refchan') && ~isempty(cfg.refchan)
    % contruct the cross-spectral-density vector of the reference channel with all MEG channels
    tmpindx = match_str(freq.labelcmb(:,1), cfg.refchan);
    refindx = match_str(freq.labelcmb(tmpindx,2), cfg.channel);
    refindx = tmpindx(refindx);
    flipref = 0;
    if isempty(refindx)
      % first look in the second column, then in the first
      tmpindx = match_str(freq.labelcmb(:,2), cfg.refchan);
      refindx = match_str(freq.labelcmb(tmpindx,1), cfg.channel);
      refindx = tmpindx(refindx);
      flipref = 1;
    end
    if length(refindx)~=Nchans
      error('The cross-spectral-density with the reference channel is not complete');
    end
    if Ntrials==1
      Cr = freq.crsspctrm(refindx, fbin);
    else
      for trial=1:Ntrials
        Cr(trial,:) = freq.crsspctrm(trial, refindx, fbin);
      end
    end
    if flipref
      Cr = conj(Cr);
    end
    % obtain the power of the reference channel
    refindx = match_str(freq.label, cfg.refchan);
    if length(refindx)<1
      error('The reference channel was not found in powspctrm');
    elseif length(refindx)>1
      error('Multiple occurences of the reference channel found in powspctrm');
    end
    if Ntrials==1
      Pr = freq.powspctrm(refindx, fbin);
    else
      for trial=1:Ntrials
        Pr(trial) = freq.powspctrm(trial, refindx, fbin);
      end
      Pr = Pr(:);   % ensure that the first dimension contains the trials
    end
  end
  
  if strcmp(cfg.dicsfix, 'yes')
    Cr = conj(Cr);
  end
  
elseif isfield(freq, 'crsspctrm')
  % this is from JMs version
  hastime    = isfield(freq, 'time');
  hasrefchan = ~isempty(cfg.refchan);
  
  % select time-frequency window of interest
  if hastime
    freq = ft_selectdata(freq, 'foilim', cfg.frequency, 'toilim', cfg.latency);
    fbin = 1;
    tbin = 1:numel(freq.time);
  else
    freq = ft_selectdata(freq, 'foilim', cfg.frequency);
    fbin = 1;
  end
  
  % convert to square csd matrix
  % think of incorporating 'quickflag' to speed up the
  % computation from fourierspectra when single trial
  % estimates are not required...
  freq = ft_checkdata(freq, 'cmbrepresentation', 'full');
  
  [dum, sensindx] = match_str(cfg.channel, freq.label);
  powspctrmindx = sensindx;
  if isempty(strfind(freq.dimord, 'rpt'))
    Ntrials = 1;
    Cf      = freq.crsspctrm(sensindx,sensindx,:,:);
    if hasrefchan,
      refindx  = match_str(freq.label, cfg.refchan);
      Cr       = freq.crsspctrm(sensindx,refindx,:,:);
      Pr       = freq.crsspctrm(refindx,refindx,:,:);
    else
      Cr = [];
      Pr = [];
    end
  elseif ~isempty(strfind(freq.dimord, 'rpt')),
    Ntrials = length(freq.cumtapcnt);
    Cf      = freq.crsspctrm(:,sensindx,sensindx,:,:);
    if hasrefchan,
      refindx  = match_str(freq.label, cfg.refchan);
      Cr       = freq.crsspctrm(:,sensindx,refindx,:,:);
      Pr       = freq.crsspctrm(:,refindx,refindx,:,:);
    else
      Cr = [];
      Pr = [];
    end
  end
else
  fprintf('computing cross-spectrum from fourier\n');
  [dum, powspctrmindx] = match_str(cfg.channel, freq.label);
  % use the fourier spectrum to compute the complete cross spectrum matrix
  Nchans = length(powspctrmindx);
  if strcmp(freq.dimord, 'chan_freq')
    error('incompatible dimord for computing CSD matrix from fourier');
  elseif strcmp(freq.dimord, 'rpt_chan_freq')
    error('incompatible dimord for computing CSD matrix from fourier');
  elseif strcmp(freq.dimord, 'rpttap_chan_freq'),
    if quickflag,
      Ntrials = 1;
    end
    Cf = zeros(Ntrials,Nchans,Nchans);
    refindx = match_str(freq.label, cfg.refchan);
    if ~isempty(refindx)
      Cr = zeros(Ntrials,Nchans,1);
      Pr = zeros(Ntrials,1,1);
    end
    
    if quickflag,
      ntap = sum(freq.cumtapcnt);
      dat  = transpose(freq.fourierspctrm(:, powspctrmindx, fbin));
      Cf(1,:,:) = (dat * ctranspose(dat)) ./ ntap;
      if ~isempty(refindx)
        ref = transpose(freq.fourierspctrm(:, refindx, fbin));
        Cr(1,:,1) = dat * ctranspose(ref) ./ ntap;
        Pr(1,1,1) = ref * ctranspose(ref) ./ ntap;
      end
    else
      freq.cumtapcnt = freq.cumtapcnt(:)';
      for k=1:Ntrials
        tapbeg = 1 + sum([0 freq.cumtapcnt(1:(k-1))]);
        tapend =     sum([0 freq.cumtapcnt(1:(k  ))]);
        ntap = freq.cumtapcnt(k);
        dat  = transpose(freq.fourierspctrm(tapbeg:tapend, powspctrmindx, fbin));
        Cf(k,:,:) = (dat * ctranspose(dat)) ./ ntap;
        if ~isempty(refindx)
          ref = transpose(freq.fourierspctrm(tapbeg:tapend, refindx, fbin));
          Cr(k,:,1) = dat * ctranspose(ref) ./ ntap;
          Pr(k,1,1) = ref * ctranspose(ref) ./ ntap;
        end
      end
    end
  else
    error('unsupported dimord for fourier cross-spectrum computation');
  end
end

% update the configuration, so that the calling function exactly knows what was selected
if ~isempty(tbin),
  % a single latency was selected in the freq structure
  cfg.latency = freq.time;
else
  cfg.latency = [];
end
cfg.frequency = freq.freq(fbin);
cfg.channel   = freq.label(powspctrmindx);

