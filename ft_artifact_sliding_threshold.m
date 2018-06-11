function [cfg, artifact,badwin] = ft_artifact_sliding_threshold(cfg, data)

% FT_ARTIFACT_SLIDING_THRESHOLD scans for time points in which the range, 
% (min-max difference) of the signal in any within a sliding time window
% exceeds a specified threshold on any channel.
% The sliding window is defined by its width (in s) and the steps (in s)
% with which it moves.
% With the 'trial' method, artifacts are whole trials. With the 'window'
% method, artifacts are time points corresponding to individual sliding
% window points.
%
% Use as
%   [cfg, artifact] = ft_artifact_threshold(cfg)
% with the configuration options
%   cfg.dataset     = string with the filename
% or
%   cfg.headerfile  = string with the filename
%   cfg.datafile    = string with the filename
%
% Alternatively you can use it as
%   [cfg, artifact] = ft_artifact_threshold(cfg, data)
%
% In both cases the configuration should also contain
%   cfg.trl        = structure that defines the data segments of interest, see FT_DEFINETRIAL
%   cfg.continuous = 'yes' or 'no' whether the file contains continuous data
%
% The following configuration options can be specified
%   cfg.artfctdef.slidethresh.method    = 'trial' or 'window'
%   cfg.artfctdef.slidethresh.channel   = cell-array with channel labels
%   cfg.artfctdef.slidethresh.bpfilter  = 'no' or 'yes'
%   cfg.artfctdef.slidethresh.bpfreq    = [0.3 30]
%   cfg.artfctdef.slidethresh.bpfiltord = 4
%
% The detection of artifacts is done according to the following settings.
%   cfg.artfctdef.slidethresh.range     = value in uV/T, default  inf
%   cfg.artfctdef.slidethresh.stdrange  = value in stdev, default  inf
%   cfg.artfctdef.slidethresh.winsize   = value in s, default = 0.1
%   cfg.artfctdef.slidethresh.winstep   = value in s, default = 0.05
%
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
if ~isfield(cfg, 'artfctdef'),          cfg.artfctdef            = [];  end
if ~isfield(cfg.artfctdef,'slidethresh'), cfg.artfctdef.slidethresh  = [];  end
if ~isfield(cfg, 'headerformat'),       cfg.headerformat         = [];  end
if ~isfield(cfg, 'dataformat'),         cfg.dataformat           = [];  end

% copy the specific configuration for this function out of the master cfg
artfctdef = cfg.artfctdef.slidethresh;

% rename some cfg fields for backward compatibility
if isfield(artfctdef, 'sgn') && ~isfield(artfctdef, 'channel')
  artfctdef.channel = artfctdef.sgn;
  artfctdef         = rmfield(artfctdef, 'sgn');
end

% set default preprocessing parameters if necessary
if ~isfield(artfctdef, 'method'),    artfctdef.method   = 'trial';    end
if ~isfield(artfctdef, 'channel'),   artfctdef.channel   = 'all';    end
if ~isfield(artfctdef, 'bpfilter'),  artfctdef.bpfilter  = 'yes';    end
if ~isfield(artfctdef, 'bpfreq'),    artfctdef.bpfreq    = [0.3 30]; end
if ~isfield(artfctdef, 'bpfiltord'), artfctdef.bpfiltord = 4;        end

% set the default artifact detection parameters
if ~isfield(artfctdef, 'range'),    artfctdef.range     = inf;        end
if ~isfield(artfctdef, 'stdrange'), artfctdef.stdrange  = inf;        end
if ~isfield(artfctdef, 'IQrange'), artfctdef.IQrange    = inf;        end
if ~isfield(artfctdef, 'winwidth'), artfctdef.winwidth  =  .1;        end
if ~isfield(artfctdef, 'winstep'),  artfctdef.winstep   =  .05;       end

% the data is either passed into the function by the user or read from file with cfg.inputfile
hasdata = exist('data', 'var');

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
  cfg.trl = data.trialinfo;
  cfg.trl(:,3) = 0;
end

% get the remaining settings
numtrl      = size(cfg.trl,1);
channel     = ft_channelselection(artfctdef.channel, hdr.label);
channelindx = match_str(hdr.label,channel);
artifact    = [];

for trlop = 1:numtrl
  if hasdata
    dat = ft_fetch_data(data,        'header', hdr, 'begsample', cfg.trl(trlop,1), 'endsample', cfg.trl(trlop,2), 'chanindx', channelindx, 'checkboundary', strcmp(cfg.continuous, 'no'));
  else
    dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', cfg.trl(trlop,1), 'endsample', cfg.trl(trlop,2), 'chanindx', channelindx, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
  end
  dat = preproc(dat, channel, offset2time(cfg.trl(trlop,3), hdr.Fs, size(dat,2)), artfctdef);

    
  t = cfg.trl(trlop,1):cfg.trl(trlop,2);
  tmin = t(1);
  tmax = t(end);
  samplestep = round(hdr.Fs*artfctdef.winstep);
  samplewidth = round(hdr.Fs*artfctdef.winwidth);
  
  wins = tmin:samplestep:tmax-samplewidth;
  wine = wins+samplewidth;
  r = NaN(numel(channel),numel(wins));
  thresh = NaN(numel(channel),numel(wins));
      
  for i = 1:numel(wins)
      win = wins(i):wine(i);
      if ~isinf(artfctdef.stdrange)
          thresh(:,i) = artfctdef.stdrange .* std(dat(:,win),[],2);
      end
      r(:,i) = range(dat(:,win),2);
      if ~isinf(artfctdef.IQrange)
         thresh(:,i) = artfctdef.IQrange .* diff(quantile(dat(:,win),[.25 .75],2),1,2);
      end 
  end
  % shift thresh to compare to previous time window
  thresh = [NaN(numel(channel),1) thresh(:,1:end-1)];
  if ~isinf(artfctdef.range)
      thresh = artfctdef.range;
  end
  badwin = r > thresh;
  n = sum(badwin(:));
  allwins = repmat(wins,size(r,1),1);
  allwine = repmat(wine,size(r,1),1);
  if strcmp(artfctdef.method,'trial')
      if n>0
          artifact(end+1,:) = [cfg.trl(trlop,1:2)];
      end
  elseif strcmp(artfctdef.method,'window')
      artifact(end+1:end+n,:) = [allwins(badwin) allwine(badwin)];
  else
      error('Unknown method');
  end
%   figure(33);
%   imagesc(r)
end

fprintf('detected %d artifacts\n', size(artifact,1));

% remember the details that were used here
cfg.artfctdef.slidethresh          = artfctdef;
cfg.artfctdef.slidethresh.wins     = wins;
cfg.artfctdef.slidethresh.trl      = cfg.trl;         % trialdefinition prior to rejection
cfg.artfctdef.slidethresh.channel  = channel;         % exact channels used for detection
cfg.artfctdef.slidethresh.artifact = artifact;        % detected artifacts

% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance
ft_postamble previous data
