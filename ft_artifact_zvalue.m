function [cfg, artifact] = ft_artifact_zvalue(cfg, data)

% FT_ARTIFACT_ZVALUE reads the interesting segments of data from file and identifies
% artifacts by means of thresholding the z-transformed value of the preprocessed raw data.
% Depending on the preprocessing options, this method will be sensitive to EOG, muscle or
% jump artifacts.  This procedure only works on continuously recorded data.
%
% Use as
%   [cfg, artifact] = ft_artifact_zvalue(cfg)
% with the configuration options
%   cfg.dataset     = string with the filename
% or
%   cfg.headerfile  = string with the filename
%   cfg.datafile    = string with the filename
% and optionally
%   cfg.headerformat
%   cfg.dataformat
%
% Alternatively you can use it as
%   [cfg, artifact] = ft_artifact_zvalue(cfg, data)
% where the input data is a structure as obtained from FT_PREPROCESSING.
%
% The required configuration settings are:
%   cfg.trl         = structure that defines the data segments of interest. See FT_DEFINETRIAL
%   cfg.continuous  = 'yes' or 'no' whether the file contains continuous data (default   = 'yes')
% and
%   cfg.artfctdef.zvalue.channel
%   cfg.artfctdef.zvalue.cutoff
%   cfg.artfctdef.zvalue.trlpadding
%   cfg.artfctdef.zvalue.fltpadding
%   cfg.artfctdef.zvalue.artpadding
%
% If you encounter difficulties with memory usage, you can use
%   cfg.memory = 'low' or 'high', whether to be memory or computationally efficient, respectively (default = 'high')
%
% The optional configuration settings (see below) are:
%   cfg.artfctdef.zvalue.artfctpeak  = 'yes' or 'no'
%   cfg.artfctdef.zvalue.interactive = 'yes' or 'no'
%
% If you specify artfctpeak='yes', the maximum value of the artifact within its range
% will be found and saved into cfg.artfctdef.zvalue.peaks.
%
% If you specify interactive='yes', a GUI will be started and you can manually
% accept/reject detected artifacts, and/or change the threshold. To control the
% graphical interface via keyboard, use the following keys:
%
%     q                 : Stop
%
%     comma             : Step to the previous artifact trial
%     a                 : Specify artifact trial to display
%     period            : Step to the next artifact trial
%
%     x                 : Step 10 trials back
%     leftarrow         : Step to the previous trial
%     t                 : Specify trial to display
%     rightarrow        : Step to the next trial
%     c                 : Step 10 trials forward
%
%     k                 : Keep trial
%     space             : Mark complete trial as artifact
%     r                 : Mark part of trial as artifact
%
%     downarrow         : Shift the z-threshold down
%     z                 : Specify the z-threshold
%     uparrow           : Shift the z-threshold down
%
% Use also, e.g. as input to DSS option of ft_componentanalysis
% cfg.artfctdef.zvalue.artfctpeakrange=[-0.25 0.25], for example to indicate range
% around peak to include, saved into cfg.artfctdef.zvalue.dssartifact. The default is
% [0 0]. Range will respect trial boundaries (i.e. be shorter if peak is near
% beginning or end of trial). Samples between trials will be removed; thus this won't
% match .sampleinfo of the data structure.
%
% Configuration settings related to the preprocessing of the data are
%   cfg.artfctdef.zvalue.lpfilter      = 'no' or 'yes'  lowpass filter
%   cfg.artfctdef.zvalue.hpfilter      = 'no' or 'yes'  highpass filter
%   cfg.artfctdef.zvalue.bpfilter      = 'no' or 'yes'  bandpass filter
%   cfg.artfctdef.zvalue.bsfilter      = 'no' or 'yes'  bandstop filter for line noise removal
%   cfg.artfctdef.zvalue.dftfilter     = 'no' or 'yes'  line noise removal using discrete fourier transform
%   cfg.artfctdef.zvalue.medianfilter  = 'no' or 'yes'  jump preserving median filter
%   cfg.artfctdef.zvalue.lpfreq        = lowpass  frequency in Hz
%   cfg.artfctdef.zvalue.hpfreq        = highpass frequency in Hz
%   cfg.artfctdef.zvalue.bpfreq        = bandpass frequency range, specified as [low high] in Hz
%   cfg.artfctdef.zvalue.bsfreq        = bandstop frequency range, specified as [low high] in Hz
%   cfg.artfctdef.zvalue.lpfiltord     = lowpass  filter order
%   cfg.artfctdef.zvalue.hpfiltord     = highpass filter order
%   cfg.artfctdef.zvalue.bpfiltord     = bandpass filter order
%   cfg.artfctdef.zvalue.bsfiltord     = bandstop filter order
%   cfg.artfctdef.zvalue.medianfiltord = length of median filter
%   cfg.artfctdef.zvalue.lpfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
%   cfg.artfctdef.zvalue.hpfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
%   cfg.artfctdef.zvalue.bpfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
%   cfg.artfctdef.zvalue.bsfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
%   cfg.artfctdef.zvalue.detrend       = 'no' or 'yes'
%   cfg.artfctdef.zvalue.demean        = 'no' or 'yes'
%   cfg.artfctdef.zvalue.baselinewindow = [begin end] in seconds, the default is the complete trial
%   cfg.artfctdef.zvalue.hilbert       = 'no' or 'yes'
%   cfg.artfctdef.zvalue.rectify       = 'no' or 'yes'
%
% The output argument "artifact" is a Nx2 matrix comparable to the
% "trl" matrix of FT_DEFINETRIAL. The first column of which specifying the
% beginsamples of an artifact period, the second column contains the
% endsamples of the artifactperiods.
%
% See also FT_REJECTARTIFACT, FT_ARTIFACT_CLIP, FT_ARTIFACT_ECG, FT_ARTIFACT_EOG,
% FT_ARTIFACT_JUMP, FT_ARTIFACT_MUSCLE, FT_ARTIFACT_THRESHOLD, FT_ARTIFACT_ZVALUE

% Copyright (C) 2003-2011, Jan-Mathijs Schoffelen & Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble loadvar data

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set default rejection parameters
cfg.headerformat                 = ft_getopt(cfg,                  'headerformat', []);
cfg.dataformat                   = ft_getopt(cfg,                  'dataformat',   []);
cfg.memory                       = ft_getopt(cfg,                  'memory',       'high');
cfg.artfctdef                    = ft_getopt(cfg,                  'artfctdef',    []);
cfg.artfctdef.zvalue             = ft_getopt(cfg.artfctdef,        'zvalue',       []);
cfg.artfctdef.zvalue.method      = ft_getopt(cfg.artfctdef.zvalue, 'method',       'all');
cfg.artfctdef.zvalue.ntrial      = ft_getopt(cfg.artfctdef.zvalue, 'ntrial',       10);
cfg.artfctdef.zvalue.channel     = ft_getopt(cfg.artfctdef.zvalue, 'channel',      {});
cfg.artfctdef.zvalue.trlpadding  = ft_getopt(cfg.artfctdef.zvalue, 'trlpadding',   0);
cfg.artfctdef.zvalue.fltpadding  = ft_getopt(cfg.artfctdef.zvalue, 'fltpadding',   0);
cfg.artfctdef.zvalue.artpadding  = ft_getopt(cfg.artfctdef.zvalue, 'artpadding',   0);
cfg.artfctdef.zvalue.interactive = ft_getopt(cfg.artfctdef.zvalue, 'interactive',  'no');
cfg.artfctdef.zvalue.cumulative  = ft_getopt(cfg.artfctdef.zvalue, 'cumulative',   'yes');
cfg.artfctdef.zvalue.artfctpeak  = ft_getopt(cfg.artfctdef.zvalue, 'artfctpeak',   'no');
cfg.artfctdef.zvalue.artfctpeakrange  = ft_getopt(cfg.artfctdef.zvalue, 'artfctpeakrange',[0 0]);

% for backward compatibility
cfg.artfctdef        = ft_checkconfig(cfg.artfctdef,        'renamed', {'blc',      'demean'});
cfg.artfctdef        = ft_checkconfig(cfg.artfctdef,        'renamed', {'blcwindow' 'baselinewindow'});
cfg.artfctdef.zvalue = ft_checkconfig(cfg.artfctdef.zvalue, 'renamed', {'sgn',      'channel'});
cfg.artfctdef.zvalue = ft_checkconfig(cfg.artfctdef.zvalue, 'renamed', {'feedback', 'interactive'});

if isfield(cfg.artfctdef.zvalue, 'artifact')
  fprintf('zvalue artifact detection has already been done, retaining artifacts\n');
  artifact = cfg.artfctdef.zvalue.artifact;
  return
end

% set feedback
cfg.feedback = ft_getopt(cfg, 'feedback',   'text');

% clear old warnings from this stack
ft_warning('-clear')

% flag whether to compute z-value per trial or not, rationale being that if
% there are fluctuations in the variance across trials (e.g. due to
% position differences in MEG measurements) which don't have to do with the artifact per se,
% the detection is compromised (although the data quality is questionable
% when there is a lot of movement to begin with).
pertrial    = strcmp(cfg.artfctdef.zvalue.method, 'trial');
demeantrial = strcmp(cfg.artfctdef.zvalue.method, 'trialdemean');
if pertrial
  if isfield(cfg.artfctdef.zvalue, 'ntrial') && cfg.artfctdef.zvalue.ntrial>0
    pertrial = cfg.artfctdef.zvalue.ntrial;
  else
    ft_error('you should specify cfg.artfctdef.zvalue.ntrial, and it should be > 0');
  end
end

% the data can be passed as input arguments or can be read from disk
hasdata = exist('data', 'var');

if ~hasdata
  % only cfg given, read data from disk
  cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
  trl = cfg.trl;
  
else
  % check whether the value for trlpadding makes sense
  if cfg.artfctdef.zvalue.trlpadding > 0
    % negative trlpadding is allowed with in-memory data
    ft_error('you cannot use positive trlpadding with in-memory data');
  end
  % check if the input data is valid for this function
  data = ft_checkdata(data, 'datatype', 'raw', 'hassampleinfo', 'yes');
  cfg = ft_checkconfig(cfg, 'forbidden', {'dataset', 'headerfile', 'datafile'});
  hdr = ft_fetch_header(data);
  trl = data.sampleinfo;
end

% set default cfg.continuous
if ~isfield(cfg, 'continuous')
  if hdr.nTrials==1
    cfg.continuous = 'yes';
  else
    cfg.continuous = 'no';
  end
end

trlpadding    = round(cfg.artfctdef.zvalue.trlpadding*hdr.Fs);
fltpadding    = round(cfg.artfctdef.zvalue.fltpadding*hdr.Fs);
artpadding    = round(cfg.artfctdef.zvalue.artpadding*hdr.Fs);
trl(:,1)      = trl(:,1) - trlpadding;       % pad the trial with some samples, in order to detect
trl(:,2)      = trl(:,2) + trlpadding;       % artifacts at the edges of the relevant trials.
if size(trl, 2) >= 3
  trl(:,3)      = trl(:,3) - trlpadding;     % the offset can ofcourse be adjusted as well
elseif hasdata
  % reconstruct offset
  for tr=1:size(trl, 1)
    % account for 0 might not be in data.time
    t0         = interp1(data.time{tr}, 1:numel(data.time{tr}), 0, 'linear', 'extrap');
    trl(tr, 3) = -t0+1 - trlpadding;
  end
else
  % assuming that the trial starts at t=0s
  trl(:, 3) = trl(:, 1);
end
trllength     = trl(:,2) - trl(:,1) + 1;     % length of each trial
numtrl        = size(trl,1);
cfg.artfctdef.zvalue.trl = trl;              % remember where we are going to look for artifacts
cfg.artfctdef.zvalue.channel = ft_channelselection(cfg.artfctdef.zvalue.channel, hdr.label);
sgnind        = match_str(hdr.label, cfg.artfctdef.zvalue.channel);
numsgn        = length(sgnind);
thresholdsum  = strcmp(cfg.artfctdef.zvalue.cumulative, 'yes');

if numsgn<1
  ft_error('no channels selected');
end

% read the data and apply preprocessing options
sumval = zeros(numsgn, 1);
sumsqr = zeros(numsgn, 1);
numsmp = zeros(numsgn, 1);
ft_progress('init', cfg.feedback, ['searching for artifacts in ' num2str(numsgn) ' channels']);
for trlop = 1:numtrl
  ft_progress(trlop/numtrl, 'searching in trial %d from %d\n', trlop, numtrl);
  
  if strcmp(cfg.memory, 'low') % store nothing in memory
    if hasdata
      dat = ft_fetch_data(data,        'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous, 'no'), 'skipcheckdata', 1);
    else
      dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
    end
    dat = preproc(dat, cfg.artfctdef.zvalue.channel, offset2time(0, hdr.Fs, size(dat,2)), cfg.artfctdef.zvalue, fltpadding, fltpadding);
    
    if trlop==1 && ~pertrial
      sumval = zeros(size(dat,1), 1);
      sumsqr = zeros(size(dat,1), 1);
      numsmp = zeros(size(dat,1), 1);
      numsgn = size(dat,1);
    elseif trlop==1 && pertrial
      sumval = zeros(size(dat,1), numtrl);
      sumsqr = zeros(size(dat,1), numtrl);
      numsmp = zeros(size(dat,1), numtrl);
      numsgn = size(dat,1);
    end
    
    if ~pertrial
      % accumulate the sum and the sum-of-squares
      sumval = sumval + nansum(dat,2);
      sumsqr = sumsqr + nansum(dat.^2,2);
      numsmp = numsmp + size(dat,2);
    else
      % store per trial the sum and the sum-of-squares
      sumval(:,trlop) = sum(dat,2);
      sumsqr(:,trlop) = sum(dat.^2,2);
      numsmp(:,trlop) = size(dat,2);
    end
  else % store all data in memory, saves computation time
    if hasdata
      dat{trlop} = ft_fetch_data(data,        'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous, 'no'), 'skipcheckdata', 1);
    else
      dat{trlop} = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
    end
    dat{trlop} = preproc(dat{trlop}, cfg.artfctdef.zvalue.channel, offset2time(0, hdr.Fs, size(dat{trlop},2)), cfg.artfctdef.zvalue, fltpadding, fltpadding);
    
    if trlop==1 && ~pertrial
      sumval = zeros(size(dat{1},1), 1);
      sumsqr = zeros(size(dat{1},1), 1);
      numsmp = zeros(size(dat{1},1), 1);
      numsgn = size(dat{1},1);
    elseif trlop==1 && pertrial
      sumval = zeros(size(dat{1},1), numtrl);
      sumsqr = zeros(size(dat{1},1), numtrl);
      numsmp = zeros(size(dat{1},1), numtrl);
      numsgn = size(dat{1},1);
    end
    
    if ~pertrial
      % accumulate the sum and the sum-of-squares
      sumval = sumval + nansum(dat{trlop},2);
      sumsqr = sumsqr + nansum(dat{trlop}.^2,2);
      numsmp = numsmp + size(dat{trlop},2);
    else
      % store per trial the sum and the sum-of-squares
      sumval(:,trlop) = sum(dat{trlop},2);
      sumsqr(:,trlop) = sum(dat{trlop}.^2,2);
      numsmp(:,trlop) = size(dat{trlop},2);
    end
  end
end % for trlop
ft_progress('close');

if pertrial>1
  sumval = ft_preproc_smooth(sumval, pertrial)*pertrial;
  sumsqr = ft_preproc_smooth(sumsqr, pertrial)*pertrial;
  numsmp = ft_preproc_smooth(numsmp, pertrial)*pertrial;
end

% compute the average and the standard deviation
datavg = sumval./numsmp;
datstd = sqrt(sumsqr./numsmp - (sumval./numsmp).^2);

if strcmp(cfg.memory, 'low')
  fprintf('\n');
end

zmax = cell(1, numtrl);
zsum = cell(1, numtrl);
zindx = cell(1, numtrl);

% create a vector that indexes the trials, or is all 1, in order
% to a per trial z-scoring, or use a static std and mean (used in lines 317
% and 328)
if pertrial
  indvec = 1:numtrl;
else
  indvec = ones(1,numtrl);
end
for trlop = 1:numtrl
  if strcmp(cfg.memory, 'low') % store nothing in memory (note that we need to preproc AGAIN... *yawn*
    fprintf('.');
    if hasdata
      dat = ft_fetch_data(data,        'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous, 'no'));
    else
      dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
    end
    dat = preproc(dat, cfg.artfctdef.zvalue.channel, offset2time(0, hdr.Fs, size(dat,2)), cfg.artfctdef.zvalue, fltpadding, fltpadding);
    zmax{trlop}  = -inf + zeros(1,size(dat,2));
    zsum{trlop}  = zeros(1,size(dat,2));
    zindx{trlop} = zeros(1,size(dat,2));
    
    nsmp          = size(dat,2);
    zdata         = (dat - datavg(:,indvec(trlop)*ones(1,nsmp)))./datstd(:,indvec(trlop)*ones(1,nsmp));  % convert the filtered data to z-values
    zsum{trlop}   = nansum(zdata,1);                   % accumulate the z-values over channels
    [zmax{trlop},ind] = max(zdata,[],1);            % find the maximum z-value and remember it
    zindx{trlop}      = sgnind(ind);                % also remember the channel number that has the largest z-value
  else
    % initialize some matrices
    zmax{trlop}  = -inf + zeros(1,size(dat{trlop},2));
    zsum{trlop}  = zeros(1,size(dat{trlop},2));
    zindx{trlop} = zeros(1,size(dat{trlop},2));
    
    nsmp          = size(dat{trlop},2);
    zdata         = (dat{trlop} - datavg(:,indvec(trlop)*ones(1,nsmp)))./datstd(:,indvec(trlop)*ones(1,nsmp));  % convert the filtered data to z-values
    zsum{trlop}   = nansum(zdata,1);                   % accumulate the z-values over channels
    [zmax{trlop},ind] = max(zdata,[],1);            % find the maximum z-value and remember it
    zindx{trlop}      = sgnind(ind);                % also remember the channel number that has the largest z-value
  end
  % This alternative code does the same, but it is much slower
  %   for i=1:size(zmax{trlop},2)
  %       if zdata{trlop}(i)>zmax{trlop}(i)
  %         % update the maximum value and channel index
  %         zmax{trlop}(i)  = zdata{trlop}(i);
  %         zindx{trlop}(i) = sgnind(sgnlop);
  %       end
  %     end
end % for trlop

if demeantrial
  for trlop = 1:numtrl
    zmax{trlop} = zmax{trlop}-mean(zmax{trlop},2);
    zsum{trlop} = zsum{trlop}-mean(zsum{trlop},2);
  end
end
%for sgnlop=1:numsgn
%  % read the data and apply preprocessing options
%  sumval = 0;
%  sumsqr = 0;
%  numsmp = 0;
%  fprintf('searching channel %s ', cfg.artfctdef.zvalue.channel{sgnlop});
%  for trlop = 1:numtrl
%    fprintf('.');
%    if hasdata
%      dat{trlop} = ft_fetch_data(data,        'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind(sgnlop), 'checkboundary', strcmp(cfg.continuous, 'no'));
%    else
%      dat{trlop} = read_data(cfg.datafile, 'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind(sgnlop), 'checkboundary', strcmp(cfg.continuous, 'no'));
%    end
%    dat{trlop} = preproc(dat{trlop}, cfg.artfctdef.zvalue.channel(sgnlop), hdr.Fs, cfg.artfctdef.zvalue, [], fltpadding, fltpadding);
%    % accumulate the sum and the sum-of-squares
%    sumval = sumval + sum(dat{trlop},2);
%    sumsqr = sumsqr + sum(dat{trlop}.^2,2);
%    numsmp = numsmp + size(dat{trlop},2);
%  end % for trlop
%
%  % compute the average and the standard deviation
%  datavg = sumval./numsmp;
%  datstd = sqrt(sumsqr./numsmp - (sumval./numsmp).^2);
%
%  for trlop = 1:numtrl
%    if sgnlop==1
%      % initialize some matrices
%      zdata{trlop} = zeros(size(dat{trlop}));
%      zmax{trlop}  = -inf + zeros(size(dat{trlop}));
%      zsum{trlop}  = zeros(size(dat{trlop}));
%      zindx{trlop} = zeros(size(dat{trlop}));
%    end
%    zdata{trlop}  = (dat{trlop} - datavg)./datstd;              % convert the filtered data to z-values
%    zsum{trlop}   = zsum{trlop} + zdata{trlop};                 % accumulate the z-values over channels
%    zmax{trlop}   = max(zmax{trlop}, zdata{trlop});             % find the maximum z-value and remember it
%    zindx{trlop}(zmax{trlop}==zdata{trlop}) = sgnind(sgnlop);   % also remember the channel number that has the largest z-value
%
%    % This alternative code does the same, but it is much slower
%    %   for i=1:size(zmax{trlop},2)
%    %       if zdata{trlop}(i)>zmax{trlop}(i)
%    %         % update the maximum value and channel index
%    %         zmax{trlop}(i)  = zdata{trlop}(i);
%    %         zindx{trlop}(i) = sgnind(sgnlop);
%    %       end
%    %     end
%  end
%  fprintf('\n');
%end % for sgnlop

for trlop = 1:numtrl
  zsum{trlop} = zsum{trlop} ./ sqrt(numsgn);
end

% always create figure
% keypress to enable keyboard uicontrol
h = figure('KeyPressFcn', @keyboard_cb);
set(h, 'visible', 'off');

opt.artcfg       = cfg.artfctdef.zvalue;
opt.artval       = {};
opt.artpadding   = artpadding;
opt.cfg          = cfg;
opt.channel      = 'artifact';
opt.hdr          = hdr;
opt.numtrl       = size(trl,1);
opt.quit         = 0;
opt.threshold    = cfg.artfctdef.zvalue.cutoff;
opt.thresholdsum = thresholdsum;
opt.trialok      = true(1,opt.numtrl); % OK by means of objective criterion
opt.keep         = zeros(1,opt.numtrl); % OK overruled by user +1 to keep, -1 to reject, start all zeros for callback to work
opt.trl          = trl;
opt.trlop        = 1;
opt.updatethreshold = true;
opt.zmax         = zmax;
opt.zsum         = zsum;

if ~thresholdsum
  opt.zval = zmax;
else
  opt.zval = zsum;
end
opt.zindx = zindx;
if ~hasdata
  opt.data = {};
else
  opt.data = data;
end

if strcmp(cfg.artfctdef.zvalue.interactive, 'yes')
  set(h, 'visible', 'on');
  set(h, 'CloseRequestFcn', @cleanup_cb);
  % give graphical feedback and allow the user to modify the threshold
  set(h, 'position', [100 200 900 400]);
  h1 = axes('position', [0.05 0.15 0.4 0.8]);
  h2 = axes('position', [0.5  0.57  0.45 0.38]);
  h3 = axes('position', [0.5  0.15  0.45 0.32]);
  opt.h1           = h1;
  opt.h2           = h2;
  opt.h3           = h3;
  
  setappdata(h, 'opt', opt);
  artval_cb(h);
  redraw_cb(h);
  
  % make the user interface elements for the data view, the order of the elements
  % here is from left to right and should match the order in the documentation
  uicontrol('tag', 'width1', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'stop',    'userdata', 'q');
  
  uicontrol('tag', 'width2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<',        'userdata', 'comma');
  uicontrol('tag', 'width1', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'artifact', 'userdata', 'a');
  uicontrol('tag', 'width2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>',        'userdata', 'period');
  
  uicontrol('tag', 'width2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<<',    'userdata', 'x');
  uicontrol('tag', 'width2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<',     'userdata', 'leftarrow');
  uicontrol('tag', 'width1', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'trial', 'userdata', 't');
  uicontrol('tag', 'width2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>',     'userdata', 'rightarrow');
  uicontrol('tag', 'width2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>>',    'userdata', 'c');
  
  uicontrol('tag', 'width3', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'keep trial',  'userdata', 'k');
  uicontrol('tag', 'width3', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'reject full', 'userdata', 'space');
  uicontrol('tag', 'width3', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'reject part', 'userdata', 'r');
  
  uicontrol('tag', 'width2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<',           'userdata', 'downarrow');
  uicontrol('tag', 'width3', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'threshold',   'userdata', 'z');
  uicontrol('tag', 'width2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>',           'userdata', 'uparrow');
  
  %uicontrol('tag', 'width2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<',       'userdata', 'control+uparrow')
  %uicontrol('tag', 'width1', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'channel', 'userdata', 'c')
  %uicontrol('tag', 'width2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>',       'userdata', 'control+downarrow')
  
  ft_uilayout(h, 'tag', 'width1', 'width', 0.10, 'height', 0.05);
  ft_uilayout(h, 'tag', 'width2', 'width', 0.05, 'height', 0.05);
  ft_uilayout(h, 'tag', 'width3', 'width', 0.12, 'height', 0.05);
  
  ft_uilayout(h, 'tag', 'width1', 'style', 'pushbutton', 'callback', @keyboard_cb);
  ft_uilayout(h, 'tag', 'width2', 'style', 'pushbutton', 'callback', @keyboard_cb);
  ft_uilayout(h, 'tag', 'width3', 'style', 'pushbutton', 'callback', @keyboard_cb);
  
  ft_uilayout(h, 'tag', 'width1', 'retag', 'viewui');
  ft_uilayout(h, 'tag', 'width2', 'retag', 'viewui');
  ft_uilayout(h, 'tag', 'width3', 'retag', 'viewui');
  ft_uilayout(h, 'tag', 'viewui', 'BackgroundColor', [0.8 0.8 0.8], 'hpos', 'auto', 'vpos', 0.005);
  
  while opt.quit==0
    uiwait(h);
    opt = getappdata(h, 'opt');
  end
  
else
  % compute the artifacts given the settings in the cfg
  setappdata(h, 'opt', opt);
  artval_cb(h);
end

h   = getparent(h);
opt = getappdata(h, 'opt');

% convert to one long vector
dum = zeros(1,max(opt.trl(:,2)));
for trlop=1:opt.numtrl
  dum(opt.trl(trlop,1):opt.trl(trlop,2)) = opt.artval{trlop};
end
artval = dum;

% find the padded artifacts and put them in a Nx2 trl-like matrix
artbeg = find(diff([0 artval])== 1);
artend = find(diff([artval 0])==-1);
artifact = [artbeg(:) artend(:)];

if strcmp(cfg.artfctdef.zvalue.artfctpeak, 'yes')
  cnt=1;
  shift=opt.trl(1,1)-1;
  for tt=1:opt.numtrl
    if tt==1
      tind{tt}=find(artifact(:,2)<opt.trl(tt,2));
    else
      tind{tt}=intersect(find(artifact(:,2)<opt.trl(tt,2)),find(artifact(:,2)>opt.trl(tt-1,2)));
    end
    artbegend=[(artifact(tind{tt},1)-opt.trl(tt,1)+1) (artifact(tind{tt},2)-opt.trl(tt,1)+1)];
    for rr=1:size(artbegend,1)
      [mx,mxnd]=max(opt.zval{tt}(artbegend(rr,1):artbegend(rr,2)));
      peaks(cnt)=artifact(tind{tt}(rr),1)+mxnd-1;
      dssartifact(cnt,1)=max(peaks(cnt)+cfg.artfctdef.zvalue.artfctpeakrange(1)*hdr.Fs,opt.trl(tt,1));
      dssartifact(cnt,2)=min(peaks(cnt)+cfg.artfctdef.zvalue.artfctpeakrange(2)*hdr.Fs,opt.trl(tt,2));
      peaks(cnt)=peaks(cnt)-shift;
      dssartifact(cnt,:)=dssartifact(cnt,:)-shift;
      cnt=cnt+1;
    end
    if tt<opt.numtrl
      shift=shift+opt.trl(tt+1,1)-opt.trl(tt,2)-1;
    end
    clear artbegend
  end
  cfg.artfctdef.zvalue.peaks=peaks';
  cfg.artfctdef.zvalue.dssartifact=dssartifact;
end

% remember the artifacts that were found
cfg.artfctdef.zvalue.artifact = artifact;

% also update the threshold
cfg.artfctdef.zvalue.cutoff   = opt.threshold;

fprintf('detected %d artifacts\n', size(artifact,1));

delete(h);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance
ft_postamble previous data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function artval_cb(h, eventdata)

opt = getappdata(h, 'opt');

artval = cell(1,opt.numtrl);
for trlop=1:opt.numtrl
  if opt.thresholdsum
    % threshold the accumulated z-values
    artval{trlop} = opt.zsum{trlop}>opt.threshold;
  else
    % threshold the max z-values
    artval{trlop} = opt.zmax{trlop}>opt.threshold;
  end
  % pad the artifacts
  artbeg = find(diff([0 artval{trlop}])== 1);
  artend = find(diff([artval{trlop} 0])==-1);
  artbeg = artbeg - opt.artpadding;
  artend = artend + opt.artpadding;
  artbeg(artbeg<1) = 1;
  artend(artend>length(artval{trlop})) = length(artval{trlop});
  for artlop=1:length(artbeg)
    artval{trlop}(artbeg(artlop):artend(artlop)) = 1;
  end
  opt.trialok(trlop) = isempty(artbeg);
end

for trlop = find(opt.keep==1 & opt.trialok==0)
  % overrule the objective criterion, i.e. keep the trial when the user
  % wants to keep it
  artval{trlop}(:) = 0;
end

for trlop = find(opt.keep<0 & opt.trialok==1)
  % if the user specifies that the trial is not OK
  % reject the whole trial if there is no extra-threshold data,
  % otherwise use the artifact as found by the thresholding
  if opt.thresholdsum && opt.keep(trlop)==-1
    % threshold the accumulated z-values
    artval{trlop} = opt.zsum{trlop}>opt.threshold;
  elseif opt.keep(trlop)==-1
    % threshold the max z-values
    artval{trlop} = opt.zmax{trlop}>opt.threshold;
  elseif opt.keep(trlop)==-2
    artval{trlop}(:) = 1;
  end
  % pad the artifacts
  artbeg = find(diff([0 artval{trlop}])== 1);
  artend = find(diff([artval{trlop} 0])==-1);
  artbeg = artbeg - opt.artpadding;
  artend = artend + opt.artpadding;
  artbeg(artbeg<1) = 1;
  artend(artend>length(artval{trlop})) = length(artval{trlop});
  if ~isempty(artbeg)
    for artlop=1:length(artbeg)
      artval{trlop}(artbeg(artlop):artend(artlop)) = 1;
    end
  else
    artval{trlop}(:) = 1;
  end
end

for trlop = find(opt.keep==-2 & opt.trialok==0)
  % if the user specifies the whole trial to be rejected define the whole
  % segment to be bad
  artval{trlop}(:) = 1;
  % pad the artifacts
  artbeg = find(diff([0 artval{trlop}])== 1);
  artend = find(diff([artval{trlop} 0])==-1);
  artbeg = artbeg - opt.artpadding;
  artend = artend + opt.artpadding;
  artbeg(artbeg<1) = 1;
  artend(artend>length(artval{trlop})) = length(artval{trlop});
  if ~isempty(artbeg)
    for artlop=1:length(artbeg)
      artval{trlop}(artbeg(artlop):artend(artlop)) = 1;
    end
  else
    artval{trlop}(:) = 1;
  end
end

opt.artval = artval;
setappdata(h, 'opt', opt);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function keyboard_cb(h, eventdata)

% If a mouseclick was made, use that value. If not, determine the key that
% corresponds to the uicontrol element that was activated.

if isa(eventdata, 'matlab.ui.eventdata.ActionData') % only the case when clicked with mouse
  curKey = get(h, 'userdata');
elseif isa(eventdata, 'matlab.ui.eventdata.KeyData') % only when key was pressed
  if isempty(eventdata.Character) && any(strcmp(eventdata.Key, {'control', 'shift', 'alt', '0'}))
    % only a modifier key was pressed
    return
  end
  if isempty(eventdata.Modifier)
    curKey = eventdata.Key;
  else
    curKey = [sprintf('%s+', eventdata.Modifier{:}) eventdata.Key];
  end
elseif isfield(eventdata, 'Key')  % only when key was pressed
  curKey = eventdata.Key;
elseif isempty(eventdata) % matlab2012b returns an empty double upon a mouse click
  curKey = get(h, 'userdata');
else
  ft_error('cannot process user input, please report this on http://bugzilla.fieldtriptoolbox.org including your MATLAB version');
end

h = getparent(h); % otherwise h is empty if isa [...].ActionData
opt = getappdata(h, 'opt');

switch strtrim(curKey)
  case 'leftarrow' % change trials
    opt.trlop = max(opt.trlop - 1, 1); % should not be smaller than 1
    setappdata(h, 'opt', opt);
    redraw_cb(h, eventdata);
  case 'x'
    opt.trlop = max(opt.trlop - 10, 1); % should not be smaller than 1
    setappdata(h, 'opt', opt);
    redraw_cb(h, eventdata);
  case 'rightarrow'
    opt.trlop = min(opt.trlop + 1, opt.numtrl); % should not be larger than the number of trials
    setappdata(h, 'opt', opt);
    redraw_cb(h, eventdata);
  case 'c'
    opt.trlop = min(opt.trlop + 10, opt.numtrl); % should not be larger than the number of trials
    setappdata(h, 'opt', opt);
    redraw_cb(h, eventdata);
  case 'uparrow' % change threshold
    opt.threshold = opt.threshold+0.5;
    opt.updatethreshold = true;
    setappdata(h, 'opt', opt);
    artval_cb(h, eventdata);
    redraw_cb(h, eventdata);
    opt = getappdata(h, 'opt'); % grab the opt-structure from the handle because it has been adjusted in the callbacks
    opt.updatethreshold = false;
    setappdata(h, 'opt', opt);
  case 'downarrow'
    opt.threshold = opt.threshold-0.5;
    opt.updatethreshold = true;
    setappdata(h, 'opt', opt);
    artval_cb(h, eventdata);
    redraw_cb(h, eventdata);
    opt = getappdata(h, 'opt'); % grab the opt-structure from the handle because it has been adjusted in the callbacks
    opt.updatethreshold = false;
    setappdata(h, 'opt', opt);
  case 'period' % change artifact
    artfctindx = find(opt.trialok == 0);
    sel        = find(artfctindx>opt.trlop);
    if ~isempty(sel)
      opt.trlop = artfctindx(sel(1));
    end
    setappdata(h, 'opt', opt);
    redraw_cb(h, eventdata);
  case 'comma'
    artfctindx = find(opt.trialok == 0);
    sel        = find(artfctindx<opt.trlop);
    if ~isempty(sel)
      opt.trlop = artfctindx(sel(end));
    end
    setappdata(h, 'opt', opt);
    redraw_cb(h, eventdata);
    %   case 'control+uparrow' % change channel
    %     if strcmp(opt.channel, 'artifact')
    %       [dum, indx] = max(opt.zval);
    %       sgnind      = opt.zindx(indx);
    %     else
    %       if ~isempty(opt.data)
    %         sgnind  = match_str(opt.channel, opt.data.label);
    %         selchan = match_str(opt.artcfg.channel, opt.channel);
    %       else
    %         sgnind  = match_str(opt.channel,   opt.hdr.label);
    %         selchan = match_str(opt.artcfg.channel, opt.channel);
    %       end
    %     end
    %     numchan = numel(opt.artcfg.channel);
    %     chansel = min(selchan+1, numchan);
    %     % convert numeric array into cell-array with channel labels
    %     opt.channel = tmpchan(chansel);
    %     setappdata(h, 'opt', opt);
    %     redraw_cb(h, eventdata);
    %   case 'c' % select channel
    %     select = match_str([opt.artcfg.channel;{'artifact'}], opt.channel);
    %     opt.channel = select_channel_list([opt.artcfg.channel;{'artifact'}], select);
    %     setappdata(h, 'opt', opt);
    %     redraw_cb(h, eventdata);
    %   case 'control+downarrow'
    %     tmpchan = [opt.artcfg.channel;{'artifact'}]; % append the 'artifact' channel
    %     chansel = match_str(tmpchan, opt.channel);
    %     chansel = max(chansel-1, 1);
    %     % convert numeric array into cell-array with channel labels
    %     opt.channel = tmpchan(chansel);
    %     setappdata(h, 'opt', opt);
    %     redraw_cb(h, eventdata);
  case 'a'
    % select the artifact to display
    response = inputdlg(sprintf('artifact trial to display'), 'specify', 1, {num2str(opt.trlop)});
    if ~isempty(response)
      artfctindx = find(opt.trialok == 0);
      sel        = str2double(response);
      sel        = min(numel(artfctindx), sel);
      sel        = max(1,                 sel);
      opt.trlop  = artfctindx(sel);
      setappdata(h, 'opt', opt);
      redraw_cb(h, eventdata);
    end
  case 'q'
    setappdata(h, 'opt', opt);
    cleanup_cb(h);
  case 't'
    % select the trial to display
    response = inputdlg(sprintf('trial to display'), 'specify', 1, {num2str(opt.trlop)});
    if ~isempty(response)
      opt.trlop = str2double(response);
      opt.trlop = min(opt.trlop, opt.numtrl); % should not be larger than the number of trials
      opt.trlop = max(opt.trlop, 1); % should not be smaller than 1
      setappdata(h, 'opt', opt);
      redraw_cb(h, eventdata);
    end
  case 'z'
    % select the threshold
    response = inputdlg('z-threshold', 'specify', 1, {num2str(opt.threshold)});
    if ~isempty(response)
      opt.threshold = str2double(response);
      opt.updatethreshold = true;
      setappdata(h, 'opt', opt);
      artval_cb(h, eventdata);
      redraw_cb(h, eventdata);
      opt = getappdata(h, 'opt'); % grab the opt-structure from the handle because it has been adjusted in the callbacks
      opt.updatethreshold = false;
      setappdata(h, 'opt', opt);
    end
  case 'k'
    opt.keep(opt.trlop) = 1;
    setappdata(h, 'opt', opt);
    artval_cb(h);
    redraw_cb(h);
  case 'r'
    % only of the trial contains a partial artifact
    if opt.trialok(opt.trlop) == 0
      opt.keep(opt.trlop) = -1;
    end
    setappdata(h, 'opt', opt);
    artval_cb(h);
    redraw_cb(h);
  case 'space'
    opt.keep(opt.trlop) = -2;
    setappdata(h, 'opt', opt);
    artval_cb(h);
    redraw_cb(h);
  case 'control+control'
    % do nothing
  case 'shift+shift'
    % do nothing
  case 'alt+alt'
    % do nothing
  otherwise
    setappdata(h, 'opt', opt);
    % this should be consistent with the help of the function
    fprintf('----------------------------------------------------------------------\n');
    fprintf('     q                 : Stop\n');
    fprintf('\n');
    fprintf('     comma             : Step to the previous artifact trial\n');
    fprintf('     a                 : Specify artifact trial to display\n');
    fprintf('     period            : Step to the next artifact trial\n');
    fprintf('\n');
    fprintf('     x                 : Step 10 trials back\n');
    fprintf('     leftarrow         : Step to the previous trial\n');
    fprintf('     t                 : Specify trial to display\n');
    fprintf('     rightarrow        : Step to the next trial\n');
    fprintf('     c                 : Step 10 trials forward\n');
    fprintf('\n');
    fprintf('     k                 : Keep trial\n');
    fprintf('     space             : Mark complete trial as artifact\n');
    fprintf('     r                 : Mark part of trial as artifact\n');
    fprintf('\n');
    fprintf('     downarrow         : Shift the z-threshold down\n');
    fprintf('     z                 : Specify the z-threshold\n');
    fprintf('     uparrow           : Shift the z-threshold down\n');
    fprintf('----------------------------------------------------------------------\n');
end
clear curKey;
uiresume(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function redraw_cb(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

% make a local copy of the relevant variables
trlop     = opt.trlop;
artval    = opt.artval{trlop};
zindx     = opt.zindx{trlop};
zval      = opt.zval{trlop};
cfg       = opt.cfg;
artcfg    = opt.artcfg;
hdr       = opt.hdr;
trl       = opt.trl;
trlpadsmp = round(artcfg.trlpadding*hdr.Fs);
channel   = opt.channel;

% determine the channel with the highest z-value to be displayed
% this is default behaviour but can be overruled in the gui
if strcmp(channel, 'artifact')
  [dum, indx] = max(zval);
  sgnind      = zindx(indx);
else
  if ~isempty(opt.data)
    sgnind = match_str(channel, opt.data.label);
  else
    sgnind = match_str(channel, hdr.label);
  end
end

if ~isempty(opt.data)
  data = ft_fetch_data(opt.data, 'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous, 'no'));
else
  data = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous, 'no'));
end

% data = preproc(data, '', hdr.Fs, artcfg, [], artcfg.fltpadding, artcfg.fltpadding);

% the string us used as title and printed in the command window
str = sprintf('trial %3d of %d, channel %s', trlop, size(trl,1), hdr.label{sgnind});
fprintf('showing %s\n', str);

%-----------------------------
% plot summary in left subplot
subplot(opt.h1); hold on;

% plot as a blue line only once
if isempty(get(opt.h1, 'children'))
  for k = 1:opt.numtrl
    xval = opt.trl(k,1):opt.trl(k,2);
    if opt.thresholdsum
      yval = opt.zsum{k};
    else
      yval = opt.zmax{k};
    end
    plot(opt.h1, xval, yval, 'linestyle', '-', 'color', 'b', 'displayname', 'data');
    xlabel('samples');
    ylabel('z-value');
  end
end
h1children = get(opt.h1, 'children');

% plot trial box
boxhandle = findall(h1children, 'displayname', 'highlight');
if isempty(boxhandle)
  % draw it
  xval = trl(opt.trlop,1):trl(opt.trlop,2);
  if opt.thresholdsum
    yval = opt.zsum{opt.trlop};
  else
    yval = opt.zmax{opt.trlop};
  end
  plot(opt.h1, xval, yval, 'linestyle', '-', 'color', 'm', 'linewidth', 2, 'displayname', 'highlight');
else
  % update it
  xval = trl(opt.trlop,1):trl(opt.trlop,2);
  if opt.thresholdsum
    yval = opt.zsum{opt.trlop};
  else
    yval = opt.zmax{opt.trlop};
  end
  set(boxhandle,  'XData', xval);
  set(boxhandle,  'YData', yval);
end

% plot as red lines the suprathreshold data points
thrhandle = findall(h1children, 'displayname', 'reddata');
if isempty(thrhandle)
  % they have to be drawn
  for k = 1:opt.numtrl
    xval = trl(k,1):trl(k,2);
    if opt.thresholdsum
      yval = opt.zsum{k};
    else
      yval = opt.zmax{k};
    end
    dum = yval<=opt.threshold;
    yval(dum) = nan;
    plot(opt.h1, xval, yval, 'linestyle', '-', 'color', [1 0 0], 'displayname', 'reddata');
  end
  hline(opt.threshold, 'color', 'r', 'linestyle', ':', 'displayname', 'threshline');
elseif ~isempty(thrhandle) && opt.updatethreshold
  % they can be updated
  for k = 1:opt.numtrl
    xval = trl(k,1):trl(k,2);
    if opt.thresholdsum
      yval = opt.zsum{k};
    else
      yval = opt.zmax{k};
    end
    dum = yval<=opt.threshold;
    yval(dum) = nan;
    set(thrhandle(k), 'XData', xval);
    set(thrhandle(k), 'YData', yval);
  end
  set(findall(h1children, 'displayname', 'threshline'), 'YData', [1 1].*opt.threshold);
end

%--------------------------------------------------
% get trial specific x-axis values and padding info
xval = ((trl(opt.trlop,1):trl(opt.trlop,2))-trl(opt.trlop,1)+trl(opt.trlop,3))./opt.hdr.Fs;
if trlpadsmp>0
  sel    = trlpadsmp:(size(data,2)-trlpadsmp);
  selpad = 1:size(data,2);
else
  sel    = 1:size(data,2);
  selpad = sel;
end

% plot data of most aberrant channel in upper subplot
subplot(opt.h2); hold on
if isempty(get(opt.h2, 'children'))
  % do the plotting
  plot(xval(selpad), data(selpad),          'color', [0.5 0.5 1], 'displayname', 'line1');
  plot(xval(sel),    data(sel),             'color', [0 0 1],     'displayname', 'line2');
  vline(xval(  1)+(trlpadsmp-1/opt.hdr.Fs), 'color', [0 0 0],     'displayname', 'vline1');
  vline(xval(end)-(trlpadsmp/opt.hdr.Fs),   'color', [0 0 0],     'displayname', 'vline2');
  data(~artval) = nan;
  plot(xval, data, 'r-', 'displayname', 'line3');
  xlabel('time(s)');
  ylabel('uV or Tesla');
  xlim([xval(1) xval(end)]);
  title(str);
else
  % update in the existing handles
  h2children = get(opt.h2, 'children');
  set(findall(h2children, 'displayname', 'vline1'), 'visible', 'off');
  set(findall(h2children, 'displayname', 'vline2'), 'visible', 'off');
  set(findall(h2children, 'displayname', 'line1'), 'XData', xval(selpad));
  set(findall(h2children, 'displayname', 'line1'), 'YData', data(selpad));
  set(findall(h2children, 'displayname', 'line2'), 'XData', xval(sel));
  set(findall(h2children, 'displayname', 'line2'), 'YData', data(sel));
  data(~artval) = nan;
  set(findall(h2children, 'displayname', 'line3'),  'XData', xval);
  set(findall(h2children, 'displayname', 'line3'),  'YData', data);
  abc2 = axis(opt.h2);
  set(findall(h2children, 'displayname', 'vline1'), 'XData', [1 1]*xval(  1)+(trlpadsmp-1/opt.hdr.Fs));
  set(findall(h2children, 'displayname', 'vline1'), 'YData', abc2(3:4));
  set(findall(h2children, 'displayname', 'vline2'), 'XData', [1 1]*xval(end)-(trlpadsmp/opt.hdr.Fs));
  set(findall(h2children, 'displayname', 'vline2'), 'YData', abc2(3:4));
  set(findall(h2children, 'displayname', 'vline1'), 'visible', 'on');
  set(findall(h2children, 'displayname', 'vline2'), 'visible', 'on');
  str = sprintf('trial %3d, channel %s', opt.trlop, hdr.label{sgnind});
  title(str);
  xlim([xval(1) xval(end)]);
end

% plot z-values in lower subplot
subplot(opt.h3); hold on;
if isempty(get(opt.h3, 'children'))
  % do the plotting
  plot(xval(selpad), zval(selpad), 'color', [0.5 0.5 1], 'displayname', 'line1b');
  plot(xval(sel),    zval(sel),    'color', [0 0 1],     'displayname', 'line2b');
  hline(opt.threshold, 'color', 'r', 'linestyle', ':', 'displayname', 'threshline');
  vline(xval(  1)+(trlpadsmp-1/opt.hdr.Fs),     'color', [0 0 0],     'displayname', 'vline1b');
  vline(xval(end)-(trlpadsmp/opt.hdr.Fs),       'color', [0 0 0],     'displayname', 'vline2b');
  zval(~artval) = nan;
  plot(xval, zval, 'r-', 'displayname', 'line3b');
  xlabel('time(s)');
  ylabel('z-value');
  xlim([xval(1) xval(end)]);
else
  % update in the existing handles
  h3children = get(opt.h3, 'children');
  set(findall(h3children, 'displayname', 'vline1b'), 'visible', 'off');
  set(findall(h3children, 'displayname', 'vline2b'), 'visible', 'off');
  set(findall(h3children, 'displayname', 'line1b'), 'XData', xval(selpad));
  set(findall(h3children, 'displayname', 'line1b'), 'YData', zval(selpad));
  set(findall(h3children, 'displayname', 'line2b'), 'XData', xval(sel));
  set(findall(h3children, 'displayname', 'line2b'), 'YData', zval(sel));
  zval(~artval) = nan;
  set(findall(h3children, 'displayname', 'line3b'),     'XData', xval);
  set(findall(h3children, 'displayname', 'line3b'),     'YData', zval);
  set(findall(h3children, 'displayname', 'threshline'), 'YData', [1 1].*opt.threshold);
  set(findall(h3children, 'displayname', 'threshline'), 'XData', xval([1 end]));
  abc = axis(opt.h3);
  set(findall(h3children, 'displayname', 'vline1b'), 'XData', [1 1]*xval(  1)+(trlpadsmp-1/opt.hdr.Fs));
  set(findall(h3children, 'displayname', 'vline1b'), 'YData', abc(3:4));
  set(findall(h3children, 'displayname', 'vline2b'), 'XData', [1 1]*xval(end)-(trlpadsmp/opt.hdr.Fs));
  set(findall(h3children, 'displayname', 'vline2b'), 'YData', abc(3:4));
  set(findall(h3children, 'displayname', 'vline1b'), 'visible', 'on');
  set(findall(h3children, 'displayname', 'vline2b'), 'visible', 'on');
  xlim([xval(1) xval(end)]);
end

setappdata(h, 'opt', opt);
uiresume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanup_cb(h, eventdata)
opt = getappdata(h, 'opt');
opt.quit = true;
setappdata(h, 'opt', opt);
uiresume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end

