function [cfg, artifact] = ft_artifact_nan(cfg, data)

% FT_ARTIFACT_NAN identifies artifacts that are indicated in the data as NaN (not a
% number) values.
%
% Use as
%   [cfg, artifact] = ft_artifact_nan(cfg, data)
% where the input data is a structure as obtained from FT_REJECTARTIFACT with the
% option cfg.artfctdef.reject='nan', or from FT_REJECTVISUAL with cfg.keeptrial='nan'
% or cfg.keepchannel='nan'.
%
% The configuration can contain
%   cfg.artfctdef.nan.channel = Nx1 cell-array with selection of channels, see FT_CHANNELSELECTION for details
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

% Copyright (C) 2017, Robert Oostenveld
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
ft_preamble debug
ft_preamble loadvar    data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set the default options
cfg.feedback              = ft_getopt(cfg, 'feedback', 'text');
cfg.artfctdef             = ft_getopt(cfg, 'artfctdef');
cfg.artfctdef.nan         = ft_getopt(cfg.artfctdef, 'nan');
cfg.artfctdef.nan.channel = ft_getopt(cfg.artfctdef.nan, 'channel', {'all'});

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
if ~isfield(cfg, 'continuous')
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
artfctdef     = cfg.artfctdef.nan;
artfctdef.trl = trl;
numtrl        = size(trl,1);
channel       = ft_channelselection(artfctdef.channel, hdr.label);
chanindx      = match_str(hdr.label, channel);
nchan         = length(chanindx);
artifact      = zeros(0,2);

ft_progress('init', cfg.feedback, ['searching for artifacts in ' num2str(nchan) ' channels']);
for trlop=1:numtrl
  ft_progress(trlop/numtrl, 'searching in trial %d from %d\n', trlop, numtrl);
  if hasdata
    dat = ft_fetch_data(data,        'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', chanindx, 'checkboundary', strcmp(cfg.continuous, 'no'));
  else
    dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', chanindx, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
  end
  
  tmp = any(isnan(dat),1);
  if any(tmp)
    % there can be multiple segments with nans
    begsample = find(diff([0 tmp])>0);
    endsample = find(diff([tmp 0])<0);
    for j=1:numel(begsample)
      artifact(end+1,:) = [begsample(j) endsample(j)] + data.sampleinfo(trlop,1) - 1;
    end
  end
end % for trlop
ft_progress('close');

% remember the details that were used here
cfg.artfctdef.nan          = artfctdef;
cfg.artfctdef.nan.artifact = artifact;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble savevar
