function [cfg, artifact] = ft_artifact_jump(cfg, data)

% FT_ARTIFACT_JUMP reads the data segments of interest from file and identifies
% SQUID jump artifacts.
%
% Use as
%   [cfg, artifact] = ft_artifact_jump(cfg)
% with the configuration options
%   cfg.dataset     = string with the filename
% or
%   cfg.headerfile  = string with the filename
%   cfg.datafile    = string with the filename
%
% Alternatively you can use it as
%   [cfg, artifact] = ft_artifact_jump(cfg, data)
%
% In both cases the configuration should also contain
%   cfg.trl        = structure that defines the data segments of interest. See FT_DEFINETRIAL
%   cfg.continuous = 'yes' or 'no' whether the file contains continuous data
%
% The data is preprocessed (again) with the following configuration parameters,
% which are optimal for identifying jump artifacts.
%   cfg.artfctdef.jump.medianfilter  = 'yes'
%   cfg.artfctdef.jump.medianfiltord = 9
%   cfg.artfctdef.jump.absdiff       = 'yes'
%
% Artifacts are identified by means of thresholding the z-transformed value
% of the preprocessed data.
%   cfg.artfctdef.jump.channel       = Nx1 cell-array with selection of channels, see FT_CHANNELSELECTION for details
%   cfg.artfctdef.jump.cutoff        = z-value at which to threshold (default = 20)
%   cfg.artfctdef.jump.trlpadding    = automatically determined based on the filter padding (cfg.padding)
%   cfg.artfctdef.jump.artpadding    = automatically determined based on the filter padding (cfg.padding)
%
% The output argument "artifact" is a Nx2 matrix comparable to the
% "trl" matrix of FT_DEFINETRIAL. The first column of which specifying the
% beginsamples of an artifact period, the second column contains the
% endsamples of the artifactperiods.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_REJECTARTIFACT, FT_ARTIFACT_CLIP, FT_ARTIFACT_ECG, FT_ARTIFACT_EOG,
% FT_ARTIFACT_JUMP, FT_ARTIFACT_MUSCLE, FT_ARTIFACT_THRESHOLD, FT_ARTIFACT_ZVALUE

% Undocumented local options:
% cfg.method

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
% ft_preamble provenance is not needed because just a call to ft_artifact_zvalue
% ft_preamble loadvar data is not needed because ft_artifact_zvalue will do this

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});

% set default rejection parameters
cfg.artfctdef               = ft_getopt(cfg,                  'artfctdef', []);
cfg.artfctdef.jump        = ft_getopt(cfg.artfctdef,        'jump',    []);
cfg.artfctdef.jump.method = ft_getopt(cfg.artfctdef.jump, 'method',    'zvalue');

if isfield(cfg.artfctdef.jump, 'artifact')
  fprintf('jump artifact detection has already been done, retaining artifacts\n');
  artifact = cfg.artfctdef.jump.artifact;
  return
end

if ~strcmp(cfg.artfctdef.jump.method, 'zvalue')
  error('jump artifact detection only works with cfg.method=''zvalue''');
end

% for backward compatibility
cfg.artfctdef.jump = ft_checkconfig(cfg.artfctdef.jump, 'renamed', {'sgn',     'channel'});
cfg.artfctdef.jump = ft_checkconfig(cfg.artfctdef.jump, 'renamed', {'passbnd', 'bpfreq'});
cfg.artfctdef.jump = ft_checkconfig(cfg.artfctdef.jump, 'renamed', {'padding', 'trlpadding'});
artpadding_oldstyle  = max(ft_getopt(cfg.artfctdef.jump, 'pretim', 0), ft_getopt(cfg.artfctdef.jump, 'psttim', 0));

% settings for preprocessing
cfg.artfctdef.jump.medianfilter  = ft_getopt(cfg.artfctdef.jump, 'medianfilter', 'yes');
cfg.artfctdef.jump.medianfiltord = ft_getopt(cfg.artfctdef.jump, 'medianfiltord', 9);
cfg.artfctdef.jump.absdiff       = ft_getopt(cfg.artfctdef.jump, 'absdiff', 'yes'); % compute abs(diff(data)), whereas the order of rectify=yes in combination with derivative=yes would be diff(abs(data)) due to the ordering in preproc

% settings for the zvalue subfunction
cfg.padding = ft_getopt(cfg, 'padding', 0);
cfg.artfctdef.jump.cutoff     = ft_getopt(cfg.artfctdef.jump, 'cutoff',     20);
cfg.artfctdef.jump.channel    = ft_getopt(cfg.artfctdef.jump, 'channel',    'MEG');
cfg.artfctdef.jump.cumulative = ft_getopt(cfg.artfctdef.jump, 'cumulative', 'no');
cfg.artfctdef.jump.trlpadding = ft_getopt(cfg.artfctdef.jump, 'trlpadding', 0.5*cfg.padding); % account for bleeding of dftfilter ringing
cfg.artfctdef.jump.fltpadding = ft_getopt(cfg.artfctdef.jump, 'fltpadding', 0.5*cfg.padding);
cfg.artfctdef.jump.artpadding = ft_getopt(cfg.artfctdef.jump, 'artpadding', 0.5*cfg.padding);

% construct a temporary configuration that can be passed onto artifact_zvalue
tmpcfg                  = cfg;
tmpcfg.artfctdef.zvalue = cfg.artfctdef.jump;
tmpcfg.artfctdef        = rmfield(tmpcfg.artfctdef, 'jump');

% call the zvalue artifact detection function, where the data is either passed
% into the function by the user or read from file with cfg.inputfile
hasdata = exist('data', 'var');
if ~hasdata
  tmpcfg = ft_checkconfig(tmpcfg, 'dataset2files', 'yes');
  tmpcfg = ft_checkconfig(tmpcfg, 'required', {'headerfile', 'datafile'});
  [tmpcfg, artifact] = ft_artifact_zvalue(tmpcfg);
else
  tmpcfg.artfctdef.zvalue.fltpadding = 0;
  warning('trlpadding and fltpadding are set to zero to avoid filter problems with NaN, see bug3193 for details');
  [tmpcfg, artifact] = ft_artifact_zvalue(tmpcfg, data);
end
cfg.artfctdef.jump = tmpcfg.artfctdef.zvalue;
