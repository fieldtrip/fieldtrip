function [cfg, artifact] = ft_artifact_muscle(cfg,data)

% FT_ARTIFACT_MUSCLE reads the data segments of interest from file and
% identifies muscle artifacts.
%
% Use as
%   [cfg, artifact] = ft_artifact_muscle(cfg)
%   required configuration options: 
%   cfg.dataset or both cfg.headerfile and cfg.datafile
% or
%   [cfg, artifact] = ft_artifact_muscle(cfg, data)
%   forbidden configuration options: 
%   cfg.dataset, cfg.headerfile and cfg.datafile
%
% In both cases the configuration should also contain:
%   cfg.trl        = structure that defines the data segments of interest. See FT_DEFINETRIAL
%   cfg.continuous = 'yes' or 'no' whether the file contains continuous data
%
% The data is preprocessed (again) with the following configuration parameters,
% which are optimal for identifying muscle artifacts:
%   cfg.artfctdef.muscle.bpfilter    = 'yes'
%   cfg.artfctdef.muscle.bpfreq      = [110 140]
%   cfg.artfctdef.muscle.bpfiltord   = 10
%   cfg.artfctdef.muscle.bpfilttype  = 'but'
%   cfg.artfctdef.muscle.hilbert     = 'yes'
%   cfg.artfctdef.muscle.boxcar      = 0.2
%
% Artifacts are identified by means of thresholding the z-transformed value
% of the preprocessed data.
%   cfg.artfctdef.muscle.channel     = Nx1 cell-array with selection of channels, see FT_CHANNELSELECTION for details
%   cfg.artfctdef.muscle.cutoff      = 4       z-value at which to threshold
%   cfg.artfctdef.muscle.trlpadding  = 0.1
%   cfg.artfctdef.muscle.fltpadding  = 0.1
%   cfg.artfctdef.muscle.artpadding  = 0.1
%
% The output argument "artifact" is a Nx2 matrix comparable to the
% "trl" matrix of FT_DEFINETRIAL. The first column of which specifying the
% beginsamples of an artifact period, the second column contains the
% endsamples of the artifactperiods.
%
% See also FT_ARTIFACT_ZVALUE, FT_REJECTARTIFACT

% Undocumented local options:
% cfg.method

% Copyright (c) 2003-2006, Jan-Mathijs Schoffelen & Robert Oostenveld
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
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
cfg = checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});

% set default rejection parameters
if ~isfield(cfg,'artfctdef'),                     cfg.artfctdef                    = [];        end
if ~isfield(cfg.artfctdef,'muscle'),              cfg.artfctdef.muscle             = [];        end
if ~isfield(cfg.artfctdef.muscle,'method'),       cfg.artfctdef.muscle.method      = 'zvalue';  end

% for backward compatibility
if isfield(cfg.artfctdef.muscle,'sgn')
  cfg.artfctdef.muscle.channel = cfg.artfctdef.muscle.sgn;
  cfg.artfctdef.muscle         = rmfield(cfg.artfctdef.muscle, 'sgn');
end

if isfield(cfg.artfctdef.muscle, 'artifact')
  fprintf('muscle artifact detection has already been done, retaining artifacts\n');
  artifact = cfg.artfctdef.muscle.artifact;
  return
end

if strcmp(cfg.artfctdef.muscle.method, 'zvalue')
  % the following settings should be supported for backward compatibility
  if isfield(cfg.artfctdef.muscle,'pssbnd'),
    cfg.artfctdef.muscle.bpfreq   = cfg.artfctdef.muscle.pssbnd;
    cfg.artfctdef.muscle.bpfilter = 'yes';
    cfg.artfctdef.muscle = rmfield(cfg.artfctdef.muscle,'pssbnd');
  end;
  dum = 0;
  if isfield(cfg.artfctdef.muscle,'pretim'),
    dum = max(dum, cfg.artfctdef.muscle.pretim);
    cfg.artfctdef.muscle = rmfield(cfg.artfctdef.muscle,'pretim');
  end
  if isfield(cfg.artfctdef.muscle,'psttim'),
    dum = max(dum, cfg.artfctdef.muscle.psttim);
    cfg.artfctdef.muscle = rmfield(cfg.artfctdef.muscle,'psttim');
  end
  if dum
    cfg.artfctdef.muscle.artpadding = max(dum);
  end
  if isfield(cfg.artfctdef.muscle,'padding'),
    cfg.artfctdef.muscle.trlpadding   = cfg.artfctdef.muscle.padding;
    cfg.artfctdef.muscle = rmfield(cfg.artfctdef.muscle,'padding');
  end

  % settings for preprocessing
  if ~isfield(cfg.artfctdef.muscle,'bpfilter'),   cfg.artfctdef.muscle.bpfilter    = 'yes';     end
  if ~isfield(cfg.artfctdef.muscle,'bpfreq'),     cfg.artfctdef.muscle.bpfreq      = [110 140]; end
  if ~isfield(cfg.artfctdef.muscle,'bpfiltord'),  cfg.artfctdef.muscle.bpfiltord   = 10;        end
  if ~isfield(cfg.artfctdef.muscle,'bpfilttype'), cfg.artfctdef.muscle.bpfilttype  = 'but';     end
  if ~isfield(cfg.artfctdef.muscle,'hilbert'),    cfg.artfctdef.muscle.hilbert     = 'yes';     end
  if ~isfield(cfg.artfctdef.muscle,'boxcar'),     cfg.artfctdef.muscle.boxcar      = 0.2;       end
  % settings for the zvalue subfunction
  if ~isfield(cfg.artfctdef.muscle,'channel'),    cfg.artfctdef.muscle.channel     = 'MEG';     end
  if ~isfield(cfg.artfctdef.muscle,'trlpadding'), cfg.artfctdef.muscle.trlpadding  = 0.1;       end
  if ~isfield(cfg.artfctdef.muscle,'fltpadding'), cfg.artfctdef.muscle.fltpadding  = 0.1;       end
  if ~isfield(cfg.artfctdef.muscle,'artpadding'), cfg.artfctdef.muscle.artpadding  = 0.1;       end
  if ~isfield(cfg.artfctdef.muscle,'cutoff'),     cfg.artfctdef.muscle.cutoff      = 4;         end
  % construct a temporary configuration that can be passed onto artifact_zvalue
  tmpcfg                  = [];
  tmpcfg.trl              = cfg.trl;
  tmpcfg.artfctdef.zvalue = cfg.artfctdef.muscle;
  if isfield(cfg, 'continuous'),   tmpcfg.continuous       = cfg.continuous;    end
  if isfield(cfg, 'dataformat'),   tmpcfg.dataformat       = cfg.dataformat;    end
  if isfield(cfg, 'headerformat'), tmpcfg.headerformat     = cfg.headerformat;  end
  % call the zvalue artifact detection function
  if nargin ==1
    cfg = checkconfig(cfg, 'dataset2files', {'yes'});
    cfg = checkconfig(cfg, 'required', {'headerfile', 'datafile'});  
    tmpcfg.datafile    = cfg.datafile;
    tmpcfg.headerfile  = cfg.headerfile;
    [tmpcfg, artifact] = ft_artifact_zvalue(tmpcfg);
  elseif nargin ==2
    cfg = checkconfig(cfg, 'forbidden', {'dataset', 'headerfile', 'datafile'});
    [tmpcfg, artifact] = ft_artifact_zvalue(tmpcfg, data);
  end
  cfg.artfctdef.muscle = tmpcfg.artfctdef.zvalue;
else
  error(sprintf('muscle artifact detection only works with cfg.method=''zvalue'''));
end

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 
