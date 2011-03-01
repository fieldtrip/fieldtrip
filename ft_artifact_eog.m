function [cfg, artifact] = ft_artifact_eog(cfg,data)

% FT_ARTIFACT_EOG reads the data segments of interest from file and
% identifies EOG artifacts.
%
% Use as
%   [cfg, artifact] = ft_artifact_eog(cfg)
%   required configuration options:
%   cfg.dataset or both cfg.headerfile and cfg.datafile
% or
%   [cfg, artifact] = ft_artifact_eog(cfg, data)
%   forbidden configuration options:
%   cfg.dataset, cfg.headerfile and cfg.datafile
%
% In both cases the configuration should also contain:
%   cfg.trl        = structure that defines the data segments of interest. See FT_DEFINETRIAL
%   cfg.continuous = 'yes' or 'no' whether the file contains continuous data
%
% The data is preprocessed (again) with the following configuration parameters,
% which are optimal for identifying EOG artifacts:
%   cfg.artfctdef.eog.bpfilter   = 'yes'
%   cfg.artfctdef.eog.bpfilttype = 'but'
%   cfg.artfctdef.eog.bpfreq     = [1 15]
%   cfg.artfctdef.eog.bpfiltord  = 4
%   cfg.artfctdef.eog.hilbert    = 'yes'
%
% Artifacts are identified by means of thresholding the z-transformed value
% of the preprocessed data.
%   cfg.artfctdef.eog.channel      = Nx1 cell-array with selection of channels, see FT_CHANNELSELECTION for details
%   cfg.artfctdef.eog.cutoff       = 4       z-value at which to threshold
%   cfg.artfctdef.eog.trlpadding   = 0.5
%   cfg.artfctdef.eog.fltpadding   = 0.1
%   cfg.artfctdef.eog.artpadding   = 0.1
%
% The output argument "artifact" is a Nx2 matrix comparable to the
% "trl" matrix of FT_DEFINETRIAL. The first column of which specifying the
% beginsamples of an artifact period, the second column contains the
% endsamples of the artifactperiods.
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following option:
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_ARTIFACT_ZVALUE, FT_REJECTARTIFACT

% Undocumented local options
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

ft_defaults

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});

% set default rejection parameters
if ~isfield(cfg,'artfctdef'),                  cfg.artfctdef                 = [];       end
if ~isfield(cfg.artfctdef,'eog'),              cfg.artfctdef.eog             = [];       end
if ~isfield(cfg.artfctdef.eog,'method'),       cfg.artfctdef.eog.method      = 'zvalue'; end
if ~isfield(cfg, 'inputfile'),                 cfg.inputfile                 = [];       end

% for backward compatibility
if isfield(cfg.artfctdef.eog,'sgn')
  cfg.artfctdef.eog.channel = cfg.artfctdef.eog.sgn;
  cfg.artfctdef.eog         = rmfield(cfg.artfctdef.eog, 'sgn');
end

if isfield(cfg.artfctdef.eog, 'artifact')
  fprintf('eog artifact detection has already been done, retaining artifacts\n');
  artifact = cfg.artfctdef.eog.artifact;
  return
end

if strcmp(cfg.artfctdef.eog.method, 'zvalue')
  % the following fields should be supported for backward compatibility
  if isfield(cfg.artfctdef.eog,'pssbnd'),
    cfg.artfctdef.eog.bpfreq   = cfg.artfctdef.eog.pssbnd;
    cfg.artfctdef.eog.bpfilter = 'yes';
    cfg.artfctdef.eog = rmfield(cfg.artfctdef.eog,'pssbnd');
  end;
  dum = 0;
  if isfield(cfg.artfctdef.eog,'pretim'),
    dum = max(dum, cfg.artfctdef.eog.pretim);
    cfg.artfctdef.eog = rmfield(cfg.artfctdef.eog,'pretim');
  end
  if isfield(cfg.artfctdef.eog,'psttim'),
    dum = max(dum, cfg.artfctdef.eog.psttim);
    cfg.artfctdef.eog = rmfield(cfg.artfctdef.eog,'psttim');
  end
  if dum
    cfg.artfctdef.eog.artpadding = max(dum);
  end
  if isfield(cfg.artfctdef.eog,'padding'),
    cfg.artfctdef.eog.trlpadding   = cfg.artfctdef.eog.padding;
    cfg.artfctdef.eog = rmfield(cfg.artfctdef.eog,'padding');
  end
  % settings for preprocessing
  if ~isfield(cfg.artfctdef.eog,'bpfilter'),   cfg.artfctdef.eog.bpfilter   = 'yes';     end
  if ~isfield(cfg.artfctdef.eog,'bpfilttype'), cfg.artfctdef.eog.bpfilttype = 'but';     end
  if ~isfield(cfg.artfctdef.eog,'bpfreq'),     cfg.artfctdef.eog.bpfreq     = [1 15];    end
  if ~isfield(cfg.artfctdef.eog,'bpfiltord'),  cfg.artfctdef.eog.bpfiltord  = 4;         end
  if ~isfield(cfg.artfctdef.eog,'hilbert'),    cfg.artfctdef.eog.hilbert    = 'yes';     end
  % settings for the zvalue subfunction
  if ~isfield(cfg.artfctdef.eog,'channel'),    cfg.artfctdef.eog.channel     = 'EOG';    end
  if ~isfield(cfg.artfctdef.eog,'trlpadding'), cfg.artfctdef.eog.trlpadding  = 0.5;      end
  if ~isfield(cfg.artfctdef.eog,'artpadding'), cfg.artfctdef.eog.artpadding  = 0.1;      end
  if ~isfield(cfg.artfctdef.eog,'fltpadding'), cfg.artfctdef.eog.fltpadding  = 0.1;      end
  if ~isfield(cfg.artfctdef.eog,'cutoff'),     cfg.artfctdef.eog.cutoff      = 4;        end
  % construct a temporary configuration that can be passed onto artifact_zvalue
  tmpcfg                  = [];
  tmpcfg.trl              = cfg.trl;
  tmpcfg.artfctdef.zvalue = cfg.artfctdef.eog;
  if isfield(cfg, 'continuous'),   tmpcfg.continuous       = cfg.continuous;    end
  if isfield(cfg, 'dataformat'),   tmpcfg.dataformat       = cfg.dataformat;    end
  if isfield(cfg, 'headerformat'), tmpcfg.headerformat     = cfg.headerformat;  end
  % call the zvalue artifact detection function
  
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
    cfg = ft_checkconfig(cfg, 'forbidden', {'dataset', 'headerfile', 'datafile'});
    [tmpcfg, artifact] = ft_artifact_zvalue(tmpcfg, data);
  else
    cfg = ft_checkconfig(cfg, 'dataset2files', {'yes'});
    cfg = ft_checkconfig(cfg, 'required', {'headerfile', 'datafile'});
    tmpcfg.datafile    = cfg.datafile;
    tmpcfg.headerfile  = cfg.headerfile;
    [tmpcfg, artifact] = ft_artifact_zvalue(tmpcfg);
  end
  cfg.artfctdef.eog  = tmpcfg.artfctdef.zvalue;
else
  error(sprintf('EOG artifact detection only works with cfg.method=''zvalue'''));
end

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');
