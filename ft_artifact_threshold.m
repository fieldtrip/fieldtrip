function [cfg, artifact] = ft_artifact_threshold(cfg,data)

% FT_ARTIFACT_THRESHOLD scans for trials in which the range, i.e. the minimum,
% the maximum or the range (min-max difference) of the signal in any
% channel exceeds a specified threshold.
%
% use as:
%   [cfg, artifact] = ft_artifact_threshold(cfg)
%   required configuration options: 
%   cfg.dataset or both cfg.headerfile and cfg.datafile
% or
%   [cfg, artifact] = ft_artifact_threshold(cfg, data)
%   forbidden configuration options: 
%   cfg.dataset, cfg.headerfile and cfg.datafile
%
% In both cases the configuration should also contain:
%   cfg.continuous                    = 'yes' or 'no' whether the file contains continuous data
%
% The following configuration options can be specified
%   cfg.artfctdef.threshold.channel   = cell-array with channel labels
%   cfg.artfctdef.threshold.bpfilter  = 'no' or 'yes'
%   cfg.artfctdef.threshold.bpfreq    = [0.3 30]
%   cfg.artfctdef.threshold.bpfiltord = 4
%
% The detection of artifacts is done according to the following settings,
% you should specify at least one of these thresholds
%   cfg.artfctdef.threshold.range     = value in uV/T, default  inf
%   cfg.artfctdef.threshold.min       = value in uV/T, default -inf
%   cfg.artfctdef.threshold.max       = value in uV/T, default  inf
%
% This function does not support partial rejections, since the whole trial
% is used to rate the minimum and maximum values. Furthermore, this
% function does not support artifact- or filterpadding.
%
% See also FT_REJECTARTIFACT

% Undocumented local options:
% cfg.inputfile
% cfg.outputfile
%
% Copyright (c) 2003, Robert Oostenveld, SMI, FCDC
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

fieldtripdefs

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'renamed',    {'ft_datatype', 'continuous'});
cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});

% set default rejection parameters for clip artifacts if necessary
if ~isfield(cfg, 'artfctdef'),          cfg.artfctdef            = [];  end
if ~isfield(cfg.artfctdef,'threshold'), cfg.artfctdef.threshold  = [];  end
if ~isfield(cfg, 'headerformat'),       cfg.headerformat         = [];  end
if ~isfield(cfg, 'dataformat'),         cfg.dataformat           = [];  end
if ~isfield(cfg, 'inputfile'),          cfg.inputfile            = [];  end
if ~isfield(cfg, 'outputfile'),         cfg.outputfile           = [];  end

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

% read the header
% depending on whether the inputfile is provided or not

hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    data = loadvar(cfg.inputfile, 'data');
    hasdata = true;
  end
end

if hasdata 
    %   isfetch = 1; 
  cfg = ft_checkconfig(cfg, 'forbidden', {'dataset', 'headerfile', 'datafile'});
  hdr = fetch_header(data);
else
    %   isfetch = 0;
  cfg = ft_checkconfig(cfg, 'dataset2files', {'yes'});
  cfg = ft_checkconfig(cfg, 'required', {'headerfile', 'datafile'});
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
end 

% set default cfg.continuous
if ~isfield(cfg, 'continuous')
    if hdr.nTrials==1
      cfg.continuous = 'yes';
    else
      cfg.continuous = 'no';
    end
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
  dat = preproc(dat, channel, hdr.Fs, artfctdef, cfg.trl(trlop,3));
  % compute the min, max and range over all channels and samples
  minval   = min(dat(:));
  maxval   = max(dat(:));
  rangeval = maxval-minval;
  % test the min, max and range against the specified thresholds
  if ~isempty(artfctdef.min) && minval<artfctdef.min
    fprintf('threshold artifact scanning: trial %d from %d exceeds min-threshold\n', trlop, numtrl);
    artifact(end+1,1:2) = cfg.trl(trlop,1:2);
  elseif ~isempty(artfctdef.max) && maxval>artfctdef.max
    fprintf('threshold artifact scanning: trial %d from %d exceeds max-threshold\n', trlop, numtrl);
    artifact(end+1,1:2) = cfg.trl(trlop,1:2);
  elseif ~isempty(artfctdef.range) && rangeval>artfctdef.range
    fprintf('threshold artifact scanning: trial %d from %d exceeds range-threshold\n', trlop, numtrl);
    artifact(end+1,1:2) = cfg.trl(trlop,1:2);
  else
    fprintf('threshold artifact scanning: trial %d from %d is ok\n', trlop, numtrl);
  end
end

% remember the details that were used here
cfg.artfctdef.threshold          = artfctdef;
cfg.artfctdef.threshold.trl      = cfg.trl;         % trialdefinition prior to rejection
cfg.artfctdef.threshold.channel  = channel;         % exact channels used for detection
cfg.artfctdef.threshold.artifact = artifact;        % detected artifacts

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id$';

if hasdata && isfield(data, 'cfg')
  % remember the configuration details of the input data
  cfg.previous = data.cfg;
end

% remember the exact configuration details in the output
data.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', data); % use the variable name "data" in the output file
end
