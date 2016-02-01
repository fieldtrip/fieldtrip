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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble init
% ft_preamble provenance is not needed because just a call to ft_artifact_zvalue
% ft_preamble loadvar data is not needed because ft_artifact_zvalue will do this

% the abort variable is set to true or false in ft_preamble_init
if abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});

% set default rejection parameters
if ~isfield(cfg,'artfctdef'),                      cfg.artfctdef                 = [];              end
if ~isfield(cfg.artfctdef,'jump'),                 cfg.artfctdef.jump            = [];              end
if ~isfield(cfg.artfctdef.jump,'method'),          cfg.artfctdef.jump.method     = 'zvalue';        end

% for backward compatibility
if isfield(cfg.artfctdef.jump,'sgn')
  cfg.artfctdef.jump.channel = cfg.artfctdef.jump.sgn;
  cfg.artfctdef.jump         = rmfield(cfg.artfctdef.jump, 'sgn');
end

if isfield(cfg.artfctdef.jump, 'artifact')
  fprintf('jump artifact detection has already been done, retaining artifacts\n');
  artifact = cfg.artfctdef.jump.artifact;
  return
end

if ~strcmp(cfg.artfctdef.jump.method, 'zvalue')
  error(sprintf('jump artifact detection only works with cfg.method=''zvalue'''));
end

% the following fields should be supported for backward compatibility
dum = 0;
if isfield(cfg.artfctdef.jump,'pretim'),
  dum = max(dum, cfg.artfctdef.jump.pretim);
  cfg.artfctdef.jump = rmfield(cfg.artfctdef.jump,'pretim');
end
if isfield(cfg.artfctdef.jump,'psttim'),
  dum = max(dum, cfg.artfctdef.jump.psttim);
  cfg.artfctdef.jump = rmfield(cfg.artfctdef.jump,'psttim');
end
if dum
  cfg.artfctdef.jump.artpadding = max(dum);
end
if isfield(cfg.artfctdef.jump,'padding'),
  cfg.artfctdef.jump.trlpadding   = cfg.artfctdef.jump.padding;
  cfg.artfctdef.jump = rmfield(cfg.artfctdef.jump,'padding');
end
% settings for preprocessing
if ~isfield(cfg.artfctdef.jump,'medianfilter'),  cfg.artfctdef.jump.medianfilter  = 'yes';        end
if ~isfield(cfg.artfctdef.jump,'medianfiltord'), cfg.artfctdef.jump.medianfiltord = 9;            end
if ~isfield(cfg.artfctdef.jump,'absdiff'),       cfg.artfctdef.jump.absdiff       = 'yes';        end  % compute abs(diff(data)), whereas the order of rectify=yes in combination with derivative=yes would be diff(abs(data)) due to the ordering in preproc
% settings for the zvalue subfunction
if ~isfield(cfg.artfctdef.jump,'cutoff'),        cfg.artfctdef.jump.cutoff     = 20;              end
if ~isfield(cfg.artfctdef.jump,'channel'),       cfg.artfctdef.jump.channel    = 'MEG';           end
if ~isfield(cfg.artfctdef.jump,'cumulative'),    cfg.artfctdef.jump.cumulative = 'no';            end
if isfield(cfg, 'padding') && cfg.padding~=0
  if ~isfield(cfg.artfctdef.jump,'trlpadding'), cfg.artfctdef.jump.trlpadding = 0.5*cfg.padding; end
  if ~isfield(cfg.artfctdef.jump,'artpadding'), cfg.artfctdef.jump.artpadding = 0.5*cfg.padding; end
  if ~isfield(cfg.artfctdef.jump,'fltpadding'), cfg.artfctdef.jump.fltpadding = 0;               end
else
  if ~isfield(cfg.artfctdef.jump,'trlpadding'), cfg.artfctdef.jump.trlpadding = 0;               end
  if ~isfield(cfg.artfctdef.jump,'artpadding'), cfg.artfctdef.jump.artpadding = 0;               end
  if ~isfield(cfg.artfctdef.jump,'fltpadding'), cfg.artfctdef.jump.fltpadding = 0;               end
end

% construct a temporary configuration that can be passed onto artifact_zvalue
tmpcfg                  = [];
tmpcfg.trl              = cfg.trl;
tmpcfg.artfctdef.zvalue = cfg.artfctdef.jump;
if isfield(cfg, 'continuous'),   tmpcfg.continuous       = cfg.continuous;    end
if isfield(cfg, 'dataformat'),   tmpcfg.dataformat       = cfg.dataformat;    end
if isfield(cfg, 'headerformat'), tmpcfg.headerformat     = cfg.headerformat;  end

% the data is either passed into the function by the user or read from file with cfg.inputfile
hasdata = exist('data', 'var');

if ~hasdata
  cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
  cfg = ft_checkconfig(cfg, 'required', {'headerfile', 'datafile'});
  tmpcfg.datafile    = cfg.datafile;
  tmpcfg.headerfile  = cfg.headerfile;
  [tmpcfg, artifact] = ft_artifact_zvalue(tmpcfg);
else
  [tmpcfg, artifact] = ft_artifact_zvalue(tmpcfg, data);
end

cfg.artfctdef.jump = tmpcfg.artfctdef.zvalue;

