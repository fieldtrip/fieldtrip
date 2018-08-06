function [cfg, artifact] = ft_artifact_threshold(cfg, data)

% FT_ARTIFACT_THRESHOLD scans for trials in which the range, i.e. the minimum, the
% maximum or the range (min-max difference) of the signal in any channel exceeds a
% specified threshold.
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
%
% The following configuration options can be specified
%   cfg.artfctdef.threshold.channel   = cell-array with channel labels
%   cfg.artfctdef.threshold.bpfilter  = 'no' or 'yes' (default = 'yes')
%   cfg.artfctdef.threshold.bpfreq    = [0.3 30]
%   cfg.artfctdef.threshold.bpfiltord = 4
%
% It is also possible to use other filter (lpfilter, hpfilter,
% bsfilter, dftfilter or medianfilter) instead of a bpfilter for
% preprocessing, see FT_PREPROCESSING
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
% Note that this function does not support artifact- or filterpadding.
%
% The output argument "artifact" is a Nx2 matrix comparable to the
% "trl" matrix of FT_DEFINETRIAL. The first column of which specifying the
% beginsamples of an artifact period, the second column contains the
% endsamples of the artifactperiods.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
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

% set default rejection parameters for clip artifacts if necessary
if ~isfield(cfg, 'artfctdef'),           cfg.artfctdef            = [];  end
if ~isfield(cfg.artfctdef, 'threshold'), cfg.artfctdef.threshold  = [];  end
if ~isfield(cfg, 'headerformat'),        cfg.headerformat         = [];  end
if ~isfield(cfg, 'dataformat'),          cfg.dataformat           = [];  end

% copy the specific configuration for this function out of the master cfg
artfctdef = cfg.artfctdef.threshold;

% rename some cfg fields for backward compatibility
if isfield(artfctdef, 'sgn') && ~isfield(artfctdef, 'channel')
  artfctdef.channel = artfctdef.sgn;
  artfctdef         = rmfield(artfctdef, 'sgn');
end
if isfield(artfctdef, 'cutoff') && ~isfield(artfctdef, 'range')
  artfctdef.range = artfctdef.cutoff;
  artfctdef       = rmfield(artfctdef, 'cutoff');
end

% set default preprocessing parameters if necessary
if ~isfield(artfctdef, 'channel'),   artfctdef.channel   = 'all';    end
if ~isfield(artfctdef, 'bpfilter'),  artfctdef.bpfilter  = 'yes';    end
if ~isfield(artfctdef, 'bpfreq'),    artfctdef.bpfreq    = [0.3 30]; end
if ~isfield(artfctdef, 'bpfiltord'), artfctdef.bpfiltord = 4;        end

% set the default artifact detection parameters
if ~isfield(artfctdef, 'range'),    artfctdef.range = inf;           end
if ~isfield(artfctdef, 'min'),      artfctdef.min =  -inf;           end
if ~isfield(artfctdef, 'max'),      artfctdef.max =   inf;           end
if ~isfield(artfctdef, 'onset'),    artfctdef.onset  = [];           end
if ~isfield(artfctdef, 'offset'),   artfctdef.offset = [];           end

% the data is either passed into the function by the user or read from file with cfg.inputfile
hasdata = exist('data', 'var');

if hasdata
  data = ft_checkdata(data, 'datatype', 'raw', 'hassampleinfo', 'yes');
end

% read the header, or get it from the input data
if ~hasdata
  cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
  cfg = ft_checkconfig(cfg, 'required', {'headerfile', 'datafile'});
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
else
  % data given as input
  data = ft_checkdata(data, 'hassampleinfo', 'yes');
  cfg = ft_checkconfig(cfg, 'forbidden', {'dataset', 'headerfile', 'datafile'});
  hdr = ft_fetch_header(data);
end

% set default cfg.continuous
if ~isfield(cfg, 'continuous')
  if hdr.nTrials==1
    cfg.continuous = 'yes';
  else
    cfg.continuous = 'no';
  end
end

if ~isfield(cfg, 'trl')
  % get it from the data itself
  cfg.trl = data.sampleinfo;
  cfg.trl(:,3) = 0;
end

% get the remaining settings
numtrl      = size(cfg.trl,1);
channel     = ft_channelselection(artfctdef.channel, hdr.label);
channelindx = match_str(hdr.label,channel);
artifact    = zeros(0,3);

if ~isempty(artfctdef.onset) || ~isempty(artfctdef.offset)
  if artfctdef.onset>0 && artfctdef.offset>0
    direction = 'up';
  elseif artfctdef.onset<0 && artfctdef.offset<0
    direction = 'down';
  else
    error('incorrect specification of onset and offset');
  end
else
  direction = 'none';
end


for trlop = 1:numtrl
  if hasdata
    dat = ft_fetch_data(data,        'header', hdr, 'begsample', cfg.trl(trlop,1), 'endsample', cfg.trl(trlop,2), 'chanindx', channelindx, 'checkboundary', strcmp(cfg.continuous, 'no'));
  else
    dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', cfg.trl(trlop,1), 'endsample', cfg.trl(trlop,2), 'chanindx', channelindx, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
  end
  
  % only do the preprocessing if there is an option that suggests to have an effect
  status = struct2cell(artfctdef);
  status = status(cellfun(@(x) ischar(x), status));
  if any(ismember(status, {'yes', 'abs', 'complex', 'real', 'imag', 'absreal', 'absimag', 'angle'}))
    dat = preproc(dat, channel, offset2time(cfg.trl(trlop,3), hdr.Fs, size(dat,2)), artfctdef);
  end
  
  % make a vector that indicates for each sample whether there is an artifact
  artval = false(1,  size(dat,2));
  artval = artval | any(dat<=artfctdef.min,1);
  artval = artval | any(dat>=artfctdef.max,1);
  
  % compute the range as the maximum of the peak-to-peak values for each channel
  ptpval = max(dat, [], 2) - min(dat, [], 2);
  if any(ptpval>=artfctdef.range)
    artval(:) = true; % mark the whole segment as bad
  end
  
  % this is when a different onset and offset are specified
  switch direction
    case 'up'
      onset  = find(diff([0 any(dat>=artfctdef.onset,1)])>0); % find the rising flank
      offset = nan(size(onset));
      for i=1:numel(onset)
        rem = dat(:,onset(i):end); % this is the remaining data following the artifact onset
        offset(i) = find(diff([any(rem<=artfctdef.offset,1) 0])>0, 1, 'first'); % find the falling flank
        offset(i) = offset(i) + onset(i);
        % add it to the other artifacts in the boolean vector
        artval(onset(i):offset(i)) = true;
      end
      
    case 'down'
      onset  = find(diff([0 any(dat<=artfctdef.onset,1)])>0); % find the rising flank
      offset = nan(size(onset));
      for i=1:numel(onset)
        rem = dat(:,onset(i):end); % this is the remaining data following the artifact onset
        offset(i) = find(diff([any(rem>=artfctdef.offset,1) 0])>0, 1, 'first'); % find the falling flank
        offset(i) = offset(i) + onset(i);
        % add it to the other artifacts in the boolean vector
        artval(onset(i):offset(i)) = true;
      end
      
    case 'none'
      % nothing to do
  end
  
  if any(artval)
    begsample = find(diff([false artval])>0) + cfg.trl(trlop,1) - 1;
    endsample = find(diff([artval false])<0) + cfg.trl(trlop,1) - 1;
    offset    = nan(size(begsample));
    
    if size(dat,1)==1
      % determine the offset of the peak value, this only works in case of a single channel
      for i=1:numel(begsample)
        seg = dat(begsample(i):endsample(i)); % get the segment of data
        if all(seg>=artfctdef.max) || strcmp(direction, 'up')
          [~, indx] = max(seg);
          offset(i) = 1 - indx; % relative to the start of the segment, 0 is the first sample, -1 is the 2nd, etc.
        elseif all(seg<=artfctdef.min) || strcmp(direction, 'down')
          [~, indx] = min(seg);
          offset(i) = 1 - indx; % relative to the start of the segment, 0 is the first sample, -1 is the 2nd, etc.
        end % if up or down
      end % for each artifact in this trial
    end % if single channel
    
    artifact  = cat(1, artifact, [begsample(:) endsample(:) offset(:)]);
  end
  
end % for trllop

if any(isnan(artifact(:,3)))
  % don't keep the offset if it cannot be determined consistently
  artifact = artifact(:,[1 2]);
end

ft_info('detected %d artifacts\n', size(artifact,1));

% remember the details that were used here
cfg.artfctdef.threshold          = artfctdef;
cfg.artfctdef.threshold.trl      = cfg.trl;         % trialdefinition prior to rejection
cfg.artfctdef.threshold.channel  = channel;         % exact channels used for detection
cfg.artfctdef.threshold.artifact = artifact;        % detected artifacts

% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance
ft_postamble previous data
