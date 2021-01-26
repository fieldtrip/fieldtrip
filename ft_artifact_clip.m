function [cfg, artifact] = ft_artifact_clip(cfg, data)

% FT_ARTIFACT_CLIP scans the data segments of interest for channels that clip, i.,e.
% channels that have a constant value for a prolonged time, often indicating that the
% signal was outside the range for the amplifier. These clipping artifacts are
% detected by the signal being completely flat for a given amount of time.
%
% Use as
%   [cfg, artifact] = ft_artifact_clip(cfg)
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
%   [cfg, artifact] = ft_artifact_clip(cfg, data)
% where the input data is a structure as obtained from FT_PREPROCESSING.
%
% In both cases the configuration should also contain
%   cfg.trl         = structure that defines the data segments of interest, see FT_DEFINETRIAL
%   cfg.continuous  = 'yes' or 'no' whether the file contains continuous data (default is automatic)
% and
%   cfg.artfctdef.clip.channel       = Nx1 cell-array with selection of channels, see FT_CHANNELSELECTION for details
%   cfg.artfctdef.clip.pretim        = pre-artifact rejection interval in seconds (default = 0)
%   cfg.artfctdef.clip.psttim        = post-artifact rejection interval in seconds (default = 0)
%   cfg.artfctdef.clip.timethreshold = number, minimum duration in seconds of a segment with consecutive identical samples to be considered as 'clipped'
%   cfg.artfctdef.clip.amplthreshold = number, minimum amplitude difference in consecutive samples to be considered as 'clipped' (default = 0)
%                                      string, percent of the amplitude range considered as 'clipped' (i.e. '1%')
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

% Copyright (C) 2005-2020, Robert Oostenveld
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
cfg = ft_checkconfig(cfg, 'renamed',    {'artfctdef.clip.thresh', 'artfctdef.clip.timethreshold'});
cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});

% set the default options
cfg.continuous      = ft_getopt(cfg, 'continuous',   []);
cfg.headerformat    = ft_getopt(cfg, 'headerformat', []);
cfg.dataformat      = ft_getopt(cfg, 'dataformat',   []);
cfg.feedback        = ft_getopt(cfg, 'feedback',   'text');
cfg.representation  = ft_getopt(cfg, 'representation', 'numeric'); % numeric or table

% set the default artifact detection parameters
cfg.artfctdef                     = ft_getopt(cfg, 'artfctdef',                    []);
cfg.artfctdef.clip                = ft_getopt(cfg.artfctdef, 'clip',               []);
cfg.artfctdef.clip.channel        = ft_getopt(cfg.artfctdef.clip, 'channel',       'all');
cfg.artfctdef.clip.timethreshold  = ft_getopt(cfg.artfctdef.clip, 'timethreshold', 0.010);
cfg.artfctdef.clip.amplthreshold  = ft_getopt(cfg.artfctdef.clip, 'amplthreshold', 0.000);
cfg.artfctdef.clip.pretim         = ft_getopt(cfg.artfctdef.clip, 'pretim',        0.000);
cfg.artfctdef.clip.psttim         = ft_getopt(cfg.artfctdef.clip, 'psttim',        0.000);

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
artfctdef     = cfg.artfctdef.clip;
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
    dat = ft_fetch_data(data,        'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', chanindx);
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
    % the threshold can be specified as a number or as a percentage
    if ischar(cfg.artfctdef.clip.amplthreshold) && endsWith(cfg.artfctdef.clip.amplthreshold, '%')
      identical = false(size(dat,1), nsample-1);
      ratio = sscanf(cfg.artfctdef.clip.amplthreshold, '%f%%');
      ratio = ratio/100;
      databsdif = abs(diff(dat(sgnlop,:), 1, 2));
      datreldif = databsdif/range(databsdif);
      identical(sgnlop,:) = datreldif <= ratio;
    else
      % detect the samples that have the same value as the previous sample
      databsdif = abs(diff(dat(sgnlop,:), 1, 2));
      identical = databsdif <= cfg.artfctdef.clip.amplthreshold;
    end
    
    % determine the number of consecutively samples for each segment where the data is identical
    count = zeros(1, nsample);
    up = find(diff([0 identical], 1, 2)== 1);
    dw = find(diff([identical 0], 1, 2)==-1) + 1;
    for k=1:length(up)
      count(up(k):dw(k)) = dw(k)-up(k)+1;
    end
    
    % detect whether there are intervals in which the number of consecutive identical samples is larger than the threshold
    thresh = (count>=artfctdef.timethreshold*hdr.Fs);
    
    begsample = find(diff([0 thresh])== 1)' + trl(trlop,1) - 1;
    endsample = find(diff([thresh 0])==-1)' + trl(trlop,1) - 1;
    channel = repmat(label(sgnlop), size(begsample));
    
    % add the pretim and psttim to the detected artifacts
    begsample = begsample - round(artfctdef.pretim * hdr.Fs);
    endsample = endsample + round(artfctdef.psttim * hdr.Fs);
    
    % remember the parts where this channel clips as artifacts
    if ~isempty(begsample)
      artifact = vertcat(artifact, table(begsample, endsample, channel));
    end
    
  end % for sgnlop
end % for trlop
ft_progress('close');

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
cfg.artfctdef.clip          = artfctdef;
cfg.artfctdef.clip.artifact = artifact;

ft_notice('detected %d artifacts\n', size(artifact,1));

% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance
ft_postamble previous data
ft_postamble savevar
