function [cfg, artifact] = ft_artifact_threshold(cfg, data)

% FT_ARTIFACT_THRESHOLD scans data segments of interest for channels in which the
% signal exceeds a specified minimum or maximum value, or in which the peak-to-peak
% range within the trial exceeds a specified threshold.
%
% Use as
%   [cfg, artifact] = ft_artifact_threshold(cfg)
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
%   [cfg, artifact] = ft_artifact_threshold(cfg, data)
% where the input data is a structure as obtained from FT_PREPROCESSING.
%
% In both cases the configuration should also contain
%   cfg.trl        = structure that defines the data segments of interest, see FT_DEFINETRIAL
%   cfg.continuous = 'yes' or 'no' whether the file contains continuous data
% and
%   cfg.artfctdef.threshold.channel   = cell-array with channel labels
%   cfg.artfctdef.threshold.bpfilter  = 'no' or 'yes' (default = 'yes')
%   cfg.artfctdef.threshold.bpfreq    = [0.3 30]
%   cfg.artfctdef.threshold.bpfiltord = 4
%
% In the same way as specifying the options for band-pass filtering, it is also
% possible to specify lpfilter, hpfilter, bsfilter, dftfilter or medianfilter, see
% FT_PREPROCESSING.
%
% The detection of artifacts is done according to the following settings,
% you should specify at least one of these thresholds
%   cfg.artfctdef.threshold.range     = value in uV or T, default  inf
%   cfg.artfctdef.threshold.min       = value in uV or T, default -inf
%   cfg.artfctdef.threshold.max       = value in uV or T, default  inf
%   cfg.artfctdef.threshold.onset     = value in uV or T, default  inf
%   cfg.artfctdef.threshold.offset    = value in uV or T, default  inf
%
% When cfg.artfctdef.threshold.range is used, the within-channel peak-to-peak range
% is checked against the specified maximum range (so not the overall range across
% channels). In this case the whole trial will be marked as an artifact.
%
% When cfg.artfctdef.threshold.onset and offset are used, the rising and falling
% flank are thresholded with different values. In case onset and offset are both
% positive, the data will be thresholded above their values. In case both onset and
% offset are negative, the data will be thresholded below their values.
%
% Note that this function does not support artifactpadding or filterpadding.
%
% The output argument "artifact" is a Nx2 matrix comparable to the "trl" matrix of
% FT_DEFINETRIAL. The first column of which specifying the beginsamples of an
% artifact period, the second column contains the endsamples of the artifactperiods.
%
% To facilitate data-handling and distributed computing, you can use
%   cfg.inputfile   =  ...
% to read the input data from a *.mat file on disk. This mat files should contain
% only a single variable named 'data', corresponding to the input structure.
%
% See also FT_REJECTARTIFACT, FT_ARTIFACT_CLIP, FT_ARTIFACT_ECG, FT_ARTIFACT_EOG,
% FT_ARTIFACT_JUMP, FT_ARTIFACT_MUSCLE, FT_ARTIFACT_THRESHOLD, FT_ARTIFACT_ZVALUE

% Copyright (C) 2003-2011, Robert Oostenveld, SMI, FCDC
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
cfg.feedback        = ft_getopt(cfg, 'feedback', 'text');
cfg.representation  = ft_getopt(cfg, 'representation', 'numeric'); % numeric or table

% set the default artifact detection parameters
cfg.artfctdef                         = ft_getopt(cfg, 'artfctdef');
cfg.artfctdef.threshold               = ft_getopt(cfg.artfctdef, 'threshold');
cfg.artfctdef.threshold.channel       = ft_getopt(cfg.artfctdef.threshold, 'channel',     'all');
cfg.artfctdef.threshold.bpfilter      = ft_getopt(cfg.artfctdef.threshold, 'bpfilter',     'yes');
cfg.artfctdef.threshold.bpfreq        = ft_getopt(cfg.artfctdef.threshold, 'bpfreq',      [0.3 30]);
cfg.artfctdef.threshold.bpfiltord     = ft_getopt(cfg.artfctdef.threshold, 'bpfiltord',   4);
cfg.artfctdef.threshold.range         = ft_getopt(cfg.artfctdef.threshold, 'range',       inf);
cfg.artfctdef.threshold.min           = ft_getopt(cfg.artfctdef.threshold, 'min',         -inf);
cfg.artfctdef.threshold.max           = ft_getopt(cfg.artfctdef.threshold, 'max',         inf);
cfg.artfctdef.threshold.onset         = ft_getopt(cfg.artfctdef.threshold, 'onset',       []);
cfg.artfctdef.threshold.offset        = ft_getopt(cfg.artfctdef.threshold, 'offset',      []);

% the data is either passed into the function by the user or read from file with cfg.inputfile
hasdata = exist('data', 'var');

% read the header, or get it from the input data
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

if ~isempty(cfg.artfctdef.threshold.onset) || ~isempty(cfg.artfctdef.threshold.offset)
  if cfg.artfctdef.threshold.onset>0 && cfg.artfctdef.threshold.offset>0
    direction = 'up';
  elseif cfg.artfctdef.threshold.onset<0 && cfg.artfctdef.threshold.offset<0
    direction = 'down';
  else
    error('incorrect specification of onset and offset');
  end
else
  direction = 'none';
end

% get the remaining settings
artfctdef     = cfg.artfctdef.threshold;
artfctdef.trl = trl;
ntrial        = size(trl,1);
label         = ft_channelselection(artfctdef.channel, hdr.label);
chanindx      = match_str(hdr.label, label);
nchan         = length(chanindx);
artifact      = table();

ft_progress('init', cfg.feedback, ['searching for artifacts in ' num2str(nchan) ' channels']);
for trlop=1:ntrial
  ft_progress(trlop/ntrial, 'searching in trial %d from %d\n', trlop, ntrial);
  if hasdata
    dat = ft_fetch_data(data,        'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', chanindx, 'checkboundary', strcmp(cfg.continuous, 'no'));
  else
    dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', chanindx, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
  end
  
  % determine the length of the data in this trial
  nsample = trl(trlop,2)-trl(trlop,1)+1;
  
  if size(trl,2)>2
    time = offset2time(trl(trlop,3), hdr.Fs, nsample);
  else
    time = offset2time(0, hdr.Fs, nsample);
  end
  
  % only do the preprocessing and filtering if there is an option that suggests to have an effect
  status = struct2cell(artfctdef);
  status = status(cellfun(@(x) ischar(x), status));
  if any(ismember(status, {'yes', 'abs', 'complex', 'real', 'imag', 'absreal', 'absimag', 'angle'}))
    dat = preproc(dat, label, time, artfctdef);
  end
  
  for sgnlop=1:nchan
    % make a vector that indicates for each sample whether it exceeds the threshold
    artval = false(1, nsample);
    artval = artval | any(dat(sgnlop,:)<=artfctdef.min,1);
    artval = artval | any(dat(sgnlop,:)>=artfctdef.max,1);
    
    % compute the range as the maximum of the peak-to-peak values for each channel
    ptpval = max(dat(sgnlop,:)) - min(dat(sgnlop,:));
    if any(ptpval>=artfctdef.range)
      artval(:) = true; % mark the whole segment as bad
    end
    
    % this is when a different onset and offset are specified
    switch direction
      case 'up'
        onset  = find(diff([0 dat(sgnlop,:)>=artfctdef.onset])>0); % find all rising flanks
        offset = nan(size(onset));
        for i=1:numel(onset)
          rem = dat(sgnlop,onset(i)+1:end); % this is the remaining data following the artifact onset
          rem = (rem<=artfctdef.offset);   % threshold for the offset
          if any(rem)
            offset(i) = find(rem, 1, 'first'); % find the falling flank
          else
            offset(i) = length(rem); % take the last sample
          end
          offset(i) = offset(i) + onset(i);
          % add it to the other artifacts in the boolean vector
          artval(onset(i):offset(i)) = true;
        end
        
      case 'down'
        onset  = find(diff([0 dat(sgnlop,:)<=artfctdef.onset])>0); % find all rising flanks
        offset = nan(size(onset));
        for i=1:numel(onset)
          rem = dat(sgnlop,onset(i)+1:end); % this is the remaining data following the artifact onset
          rem = (rem>=artfctdef.offset);
          if any(rem)
            offset(i) = find(rem, 1, 'first'); % find the falling flank
          else
            offset(i) = length(rem); % take the last sample
          end
          offset(i) = offset(i) + onset(i);
          % add it to the other artifacts in the boolean vector
          artval(onset(i):offset(i)) = true;
        end
        
      case 'none'
        % nothing to do
    end
    
    % to avoid confusion with the offset that is used further down
    clear onset offset
    
    begsample = find(diff([0 artval])>0)';
    endsample = find(diff([artval 0])<0)';
    offset    = nan(size(begsample)); % the offset of the peak relative to the segment, just like in FT_DEFINETRIAL
    channel   = repmat(label(sgnlop), size(begsample));
    
    % determine the sample at which the signal peaks
    for i=1:numel(begsample)
      seg = dat(sgnlop,begsample(i):endsample(i)); % get the segment of data
      if all(seg>=artfctdef.max) || strcmp(direction, 'up')
        [dum, indx] = max(seg);
        offset(i)   = 1 - indx; % relative to the start of the segment, 0 is the first sample, -1 is the 2nd, etc.
      elseif all(seg<=artfctdef.min) || strcmp(direction, 'down')
        [dum, indx] = min(seg);
        offset(i)   = 1 - indx; % relative to the start of the segment, 0 is the first sample, -1 is the 2nd, etc.
      end % if up or down
    end % for each artifact in this trial
    
    % express them relative to the start of the data, not the start of the trial
    begsample = begsample + trl(trlop,1) - 1;
    endsample = endsample + trl(trlop,1) - 1;
    
    % remember the parts where this channel exceeds the threshold as artifacts
    if ~isempty(begsample)
      artifact = vertcat(artifact, table(begsample, endsample, offset, channel));
    end
    
  end % for sgnlop
end % for trlop
ft_progress('close');

if strcmp(cfg.representation, 'numeric') && istable(artifact)
  if isempty(artifact)
    % an empty table does not have columns
    artifact = zeros(0,3);
  else
    % convert the table to a numeric array with the columns begsample, endsample and offset
    artifact = table2array(artifact(:,1:3));
  end
elseif strcmp(cfg.representation, 'table') && isnumeric(artifact)
  if isempty(artifact)
    % an empty table does not have columns
    artifact = table();
  else
    % convert the numeric array to a table with the columns begsample, endsample and offset
    begsample = artifact(:,1);
    endsample = artifact(:,2);
    offset    = artifact(:,3);
    artifact = table(begsample, endsample, offset);
  end
end

% remember the details that were used here and store the detected artifacts
cfg.artfctdef.threshold          = artfctdef;
cfg.artfctdef.threshold.artifact = artifact;

ft_notice('detected %d artifacts\n', size(artifact,1));

% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance
ft_postamble previous data
ft_postamble savevar
