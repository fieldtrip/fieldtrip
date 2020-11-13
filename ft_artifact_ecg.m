function [cfg, artifact] = ft_artifact_ecg(cfg, data)

% FT_ARTIFACT_ECG performs a peak-detection on the ECG-channel and identifies the
% windows around the QRS peak as artifacts. Using FT_REJECTARTIFACT you can remove
% these windows from your data, or using FT_REMOVETEMPLATEARTIFACT you can subtract
% an averaged template artifact from your data.
%
% Use as
%   [cfg, artifact] = ft_artifact_ecg(cfg)
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
%   [cfg, artifact] = ft_artifact_ecg(cfg, data)
% where the input data is a structure as obtained from FT_PREPROCESSING.
%
% In both cases the configuration should also contain
%   cfg.trl        = structure that defines the data segments of interest. See FT_DEFINETRIAL
%   cfg.continuous = 'yes' or 'no' whether the file contains continuous data
% and
%   cfg.artfctdef.ecg.channel = Nx1 cell-array with selection of channels, see FT_CHANNELSELECTION for details
%   cfg.artfctdef.ecg.pretim  = pre-artifact rejection interval in seconds (default = 0.05)
%   cfg.artfctdef.ecg.psttim  = post-artifact rejection interval in seconds (default = 0.3)
%   cfg.artfctdef.ecg.cutoff  = peak threshold (default = 3)
%   cfg.artfctdef.ecg.inspect = Nx1 list of channels which will be shown as a QRS-locked average
%
% The output argument "artifact" is a Nx2 matrix comparable to the "trl" matrix of
% FT_DEFINETRIAL. The first column of which specifying the begin samples of an
% artifact period, the second column contains the end samples of the QRS periods.
%
% To facilitate data-handling and distributed computing, you can use
%   cfg.inputfile   =  ...
% to read the input data from a *.mat file on disk. This mat files should contain
% only a single variable named 'data', corresponding to the input structure.
%
% See also FT_REJECTARTIFACT, FT_REMOVETEMPLATEARTIFACT, FT_ARTIFACT_CLIP, FT_ARTIFACT_ECG,
% FT_ARTIFACT_EOG, FT_ARTIFACT_JUMP, FT_ARTIFACT_MUSCLE, FT_ARTIFACT_THRESHOLD,
% FT_ARTIFACT_ZVALUE

% Copyright (C) 2005-2011, Jan-Mathijs Schoffelen
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

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});

% set the default options
cfg.continuous      = ft_getopt(cfg, 'continuous',   []);
cfg.headerformat    = ft_getopt(cfg, 'headerformat', []);
cfg.dataformat      = ft_getopt(cfg, 'dataformat',   []);
cfg.feedback        = ft_getopt(cfg, 'feedback',   'text');
cfg.representation  = ft_getopt(cfg, 'representation', 'numeric'); % numeric or table

% set the default artifact detection parameters
cfg.artfctdef               = ft_getopt(cfg, 'artfctdef',              []);
cfg.artfctdef.ecg           = ft_getopt(cfg.artfctdef, 'ecg',          []);
cfg.artfctdef.ecg.channel   = ft_getopt(cfg.artfctdef.ecg, 'channel',  {'ECG'});
cfg.artfctdef.ecg.method    = ft_getopt(cfg.artfctdef.ecg, 'method',   'zvalue');
cfg.artfctdef.ecg.cutoff    = ft_getopt(cfg.artfctdef.ecg, 'cutoff',   3);
cfg.artfctdef.ecg.padding   = ft_getopt(cfg.artfctdef.ecg, 'padding',  0.5);
cfg.artfctdef.ecg.inspect   = ft_getopt(cfg.artfctdef.ecg, 'inspect',  {'MLT' 'MRT'});
cfg.artfctdef.ecg.pretim    = ft_getopt(cfg.artfctdef.ecg, 'pretim',   0.05);
cfg.artfctdef.ecg.psttim    = ft_getopt(cfg.artfctdef.ecg, 'psttim',   0.3);
cfg.artfctdef.ecg.mindist   = ft_getopt(cfg.artfctdef.ecg, 'mindist',  0.5);
cfg.artfctdef.ecg.feedback  = ft_getopt(cfg.artfctdef.ecg, 'feedback', 'yes');

if ~strcmp(cfg.artfctdef.ecg.method, 'zvalue')
  ft_error('method "%s" is not applicable', cfg.artfctdef.ecg.method);
end

% the data is either passed into the function by the user or read from file with cfg.inputfile
hasdata = exist('data', 'var');

if ~hasdata
  cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
  cfg = ft_checkconfig(cfg, 'required', {'headerfile', 'datafile'});
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
else
  data = ft_checkdata(data, 'datatype', 'raw', 'hassampleinfo', 'yes');
  cfg  = ft_checkconfig(cfg, 'forbidden', {'dataset', 'headerfile', 'datafile'});
  hdr  = ft_fetch_header(data);
end

% set default cfg.continuous
if isempty(cfg.continuous)
  if hdr.nTrials==1
    cfg.continuous = 'yes';
  else
    cfg.continuous = 'no';
  end
end

% get the specification of the data segments that should be scanned for artifacts
if ~isfield(cfg, 'trl') && hasdata
  trl = data.sampleinfo;
  for k = 1:numel(data.trial)
    trl(k,3) = time2offset(data.time{k}, data.fsample);
  end
elseif isfield(cfg, 'trl') && ischar(cfg.trl)
  trl = loadvar(cfg.trl, 'trl');
elseif isfield(cfg, 'trl') && isnumeric(cfg.trl)
  trl = cfg.trl;
else
  ft_error('cannot determine which segments of data to scan for artifacts');
end

% get the remaining settings
artfctdef         = cfg.artfctdef.ecg;
artfctdef.trl     = trl;
artfctdef.demean  = 'yes';
artfctdef.channel = ft_channelselection(artfctdef.channel, hdr.label);
chanindx          = match_str(hdr.label, artfctdef.channel);
nchan             = length(chanindx);
fltpadding        = 0;
numtrl            = size(trl,1);

if length(chanindx)<1
  ft_error('no ECG channels selected');
elseif length(chanindx)>1
  ft_error('only one ECG channel can be selected');
end

% these are the settings for filtering, rectifying, etc.
fltcfg = removefields(artfctdef, {'pretim', 'psttim', 'method', 'cutoff', 'inspect'});

ft_progress('init', cfg.feedback, ['searching for artifacts in ' num2str(nchan) ' channels']);
for trlop=1:numtrl
  ft_progress(trlop/numtrl, 'searching in trial %d from %d\n', trlop, numtrl);
  if hasdata
    dat = ft_fetch_data(data,        'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', chanindx, 'checkboundary', strcmp(cfg.continuous, 'no'));
  else
    dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', chanindx, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
  end
  if size(trl,2)>2
    fltdata.time{trlop} = offset2time(trl(trlop,3), hdr.Fs, size(dat,2));
  else
    fltdata.time{trlop} = offset2time(0, hdr.Fs, size(dat,2));
  end
  fltdata.trial{trlop} = preproc(dat, artfctdef.channel, fltdata.time{trlop}, fltcfg, fltpadding, fltpadding);
  fltdata.trial{trlop} = fltdata.trial{trlop}.^2;
end % for trlop
ft_progress('close');

tmp   = cell2mat(fltdata.trial);
stmp  = std(tmp, 0, 2);
mtmp  = mean(tmp, 2);
Nsmp  = max(trl(:,2));
trace = zeros(1, Nsmp);

% standardise the ecg
for trlop=1:numtrl
  trace(trl(trlop,1):trl(trlop,2)) = (fltdata.trial{trlop}-mtmp)./stmp;
end

accept = strcmp(cfg.artfctdef.ecg.feedback, 'no');
while accept == 0
  h = figure;
  plot(trace);zoom;
  hold on;
  plot([1 Nsmp], [artfctdef.cutoff artfctdef.cutoff], 'r:');
  hold off;
  xlabel('samples');
  ylabel('zscore');
  
  fprintf(['\ncurrent %s threshold = %1.3f'], artfctdef.method, artfctdef.cutoff);
  response = input('\nkeep the current value (y/n) ?\n', 's');
  switch response
    case 'n'
      oldcutoff = artfctdef.cutoff;
      artfctdef.cutoff = input('\nenter new value \n');
    case 'y'
      oldcutoff = artfctdef.cutoff;
      accept = 1;
    otherwise
      ft_warning('unrecognised response, assuming no');
      oldcutoff = artfctdef.cutoff;
      artfctdef.cutoff = input('\nenter new value \n');
  end
  close
end

% detect peaks which are at least half a second apart and store
% the indices of the qrs-complexes in the artifact-configuration
mindist       = round(cfg.artfctdef.ecg.mindist.*hdr.Fs);
[pindx, pval] = peakdetect2(trace, artfctdef.cutoff, mindist);
%sel           = find(standardise(pval,2)<2);
%pindx         = pindx(sel);
%pval          = pval(sel);
artfctdef.qrs = pindx;

% create a trial around each QRS peak
trl = [];
trl(:,1) = pindx(:) - round(artfctdef.padding*(hdr.Fs))  ;
trl(:,2) = pindx(:) + round(artfctdef.padding*(hdr.Fs))-1;
trl(:,3) = -round(artfctdef.padding*(hdr.Fs));
trl(trl(:,1)<1,:) = [];
trl(trl(:,2)>hdr.nSamples.*hdr.nTrials,:) = [];

% ---------------------
% compute a QRS-triggered average
% FIXME, at present this only works for continuous data: the assumption can be made that all trials are equally long.
sgn    = ft_channelselection(artfctdef.inspect, hdr.label);
megind = match_str(hdr.label, sgn);
chanindx = [megind(:); chanindx];
dat    = zeros(length(chanindx), trl(1,2)-trl(1,1)+1);
numtrl   = size(trl,1);

if ~isempty(chanindx)
  ntrlok = 0;
  for trlop=1:numtrl
    if ~hasdata
      fprintf('reading and preprocessing heartbeat %d of %d\n', trlop, numtrl);
      dum = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', chanindx, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
      dat = dat + ft_preproc_baselinecorrect(dum);
      ntrlok = ntrlok + 1;
    elseif hasdata
      fprintf('preprocessing heartbeat %d of %d\n', trlop, numtrl);
      dum = ft_fetch_data(data, 'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', chanindx, 'checkboundary', strcmp(cfg.continuous, 'no'), 'docheck', 0);
      if any(~isfinite(dum(:)))
        % do not add this segment to the sum
      else
        dat    = dat + ft_preproc_baselinecorrect(dum);
        ntrlok = ntrlok + 1;
      end
    end
  end
end

dat  = dat./ntrlok;
time = offset2time(trl(1,3), hdr.Fs, size(dat,2));
tmp  = dat(1:end-1,:);
mdat = max(abs(tmp(:)));

acceptpre = strcmp(cfg.artfctdef.ecg.feedback, 'no');
acceptpst = strcmp(cfg.artfctdef.ecg.feedback, 'no');
while acceptpre == 0 || acceptpst == 0
  h = figure;
  subplot(2,1,1); plot(time, dat(end, :));
  abc = axis;
  axis([time(1) time(end) abc(3:4)]);
  subplot(2,1,2);
  axis([time(1) time(end) -1.1*mdat 1.1*mdat]);
  xpos   = -artfctdef.pretim;
  ypos   = -1.05*mdat;
  width  = artfctdef.pretim + artfctdef.psttim;
  height = 2.1*mdat;
  rectangle('Position', [xpos ypos width height], 'FaceColor', 'r');
  hold on; plot(time, dat(1:end-1, :), 'b');
  
  if acceptpre == 0
    fprintf(['\ncurrent pre-peak interval = %1.3f'], artfctdef.pretim);
    response = input('\nkeep the current value (y/n) ?\n', 's');
    switch response
      case 'n'
        oldpretim = artfctdef.pretim;
        artfctdef.pretim = input('\nenter new value \n');
      case 'y'
        oldpretim = artfctdef.pretim;
        acceptpre = 1;
      otherwise
        ft_warning('unrecognised response, assuming no');
        oldpretim = artfctdef.pretim;
    end
  end
  if acceptpst == 0 && acceptpre == 1
    fprintf(['\ncurrent post-peak interval = %1.3f'], artfctdef.psttim);
    response = input('\nkeep the current value (y/n) ?\n', 's');
    switch response
      case 'n'
        oldpsttim = artfctdef.psttim;
        artfctdef.psttim = input('\nenter new value \n');
      case 'y'
        oldpsttim = artfctdef.psttim;
        acceptpst = 1;
      otherwise
        ft_warning('unrecognised response, assuming no');
        oldpsttim = artfctdef.psttim;
    end
  end
  close
end

artifact(:,1) = trl(:,1) - trl(:,3) - round(artfctdef.pretim*hdr.Fs);
artifact(:,2) = trl(:,1) - trl(:,3) + round(artfctdef.psttim*hdr.Fs);

if strcmp(cfg.representation, 'numeric') && istable(artifact)
  if isempty(artifact)
    % an empty table does not have columns
    artifact = zeros(0,2);
  else
    % convert the table to a numeric array with the columns begsample and endsample
    artifact = table2array(artifact(:,1:2));
  end
elseif strcmp(cfg.representation, 'table') && isnumeric(artifact)
  if isempty(artifact)
    % an empty table does not have columns
    artifact = table();
  else
    % convert the numeric array to a table with the columns begsample and endsample
    begsample = artifact(:,1);
    endsample = artifact(:,2);
    artifact = table(begsample, endsample);
  end
end

% remember the details that were used here and store the detected artifacts
cfg.artfctdef.ecg          = artfctdef;
cfg.artfctdef.ecg.artifact = artifact;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance
ft_postamble previous data
ft_postamble savevar
