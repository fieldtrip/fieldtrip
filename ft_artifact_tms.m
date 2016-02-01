function [cfg, artifact] = ft_artifact_tms(cfg, data)

% FT_ARTIFACT_TMS reads the data segments of interest from file and
% identifies tms artifacts.
%
% Use as
%   [cfg, artifact] = ft_artifact_tms(cfg)
% with the configuration options
%   cfg.dataset     = string with the filename
% or
%   cfg.headerfile  = string with the filename
%   cfg.datafile    = string with the filename
%
% Alternatively you can use it as
%   [cfg, artifact] = ft_artifact_tms(cfg, data)
%
% In both cases the configuration should also contain
%   cfg.trl         = structure that defines the data segments of interest. See FT_DEFINETRIAL
%   cfg.continuous  = 'yes' or 'no' whether the file contains continuous data (default   = 'yes')
%   cfg.method      = 'detect', TMS-artifacts are detected by preprocessing
%                     the data to be sensitive to transient high gradients, typical for
%                     TMS-pulses.
%                     'marker', TMS-artifact onset and offsets are based on
%                     markers written in the EEG.
%   cfg.prestim     = scalar, time in seconds prior to onset of detected
%                     event to mark as artifactual (default = 0.005 seconds)
%   cfg.poststim    = scalar, time in seconds post onset of detected even to
%                     mark as artifactual (default = 0.010 seconds)
%
% METHOD SPECIFIC OPTIONS AND DESCRIPTIONS
%
% DETECT
% The data is preprocessed (again) with the following configuration parameters,
% which are optimal for identifying tms artifacts. This acts as a wrapper
% around ft_artifact_zvalue
%   cfg.artfctdef.tms.derivative  = 'yes'
%
% Artifacts are identified by means of thresholding the z-transformed value
% of the preprocessed data.
%   cfg.artfctdef.tms.channel     = Nx1 cell-array with selection of channels, see FT_CHANNELSELECTION for details
%   cfg.artfctdef.tms.cutoff      = z-value at which to threshold (default = 4)
%   cfg.artfctdef.tms.trlpadding  = 0.1
%   cfg.artfctdef.tms.fltpadding  = 0.1
%   cfg.artfctdef.tms.artpadding  = 0.01 (Be aware that if one artifact
%   falls within this specified range of another artifact, both artifact
%   will be counted as one. Depending on cfg.prestim and cfg.poststim you
%   may not mark enough data as artifactual.)
%
% MARKER
% This method acts as a wrapper around FT_DEFINETRIAL to determine on- and
% offsets of TMS pulses by reading markers in the EEG.
%   cfg.trialfun            = function name, see below (default = 'ft_trialfun_general')
%   cfg.trialdef.eventtype  = 'string'
%   cfg.trialdef.eventvalue = number, string or list with numbers or strings
%
% The cfg.trialfun option is a string containing the name of a function
% that you wrote yourself and that FT_ARTIFACT_TMS will call. The
% function should take the cfg-structure as input and should give a
% NxM matrix with M equal to or larger than 3) in the same format as
% "trl" as the output. You can add extra custom fields to the
% configuration structure to pass as arguments to your own trialfun.
% Furthermore, inside the trialfun you can use the FT_READ_EVENT
% function to get the event information from your data file.
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
% See also FT_REJECTARTIFACT, FT_ARTIFACT_CLIP, FT_ARTIFACT_ECG, FT_ARTIFACT_EOG,
% FT_ARTIFACT_JUMP, FT_ARTIFACT_MUSCLE, FT_ARTIFACT_THRESHOLD, FT_ARTIFACT_ZVALUE

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

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});
cfg = ft_checkconfig(cfg, 'required', 'method');
cfg = ft_checkconfig(cfg, 'allowedval',{'method','detect','marker'});

% set default rejection parameters
if ~isfield(cfg,'artfctdef'),                       cfg.artfctdef                   = [];        end
if ~isfield(cfg,'method'),                          cfg.method                      = 'detect';  end
if ~isfield(cfg.artfctdef,'tms'),                   cfg.artfctdef.tms               = [];        end
if ~isfield(cfg,'prestim'),                         cfg.prestim                     = 0.005;     end
if ~isfield(cfg,'poststim'),                        cfg.poststim                    = 0.010;     end

if isfield(cfg.artfctdef.tms, 'artifact')
  fprintf('tms artifact detection has already been done, retaining artifacts\n');
  artifact = cfg.artfctdef.tms.artifact;
  return
end

switch cfg.method
  case 'detect'
    % settings for preprocessing
    if ~isfield(cfg.artfctdef.tms,'derivative'),   cfg.artfctdef.tms.derivative   = 'yes';    end
    % settings for the zvalue subfunction
    if ~isfield(cfg.artfctdef.tms,'method'),     cfg.artfctdef.tms.method      = 'zvalue';    end
    if ~isfield(cfg.artfctdef.tms,'channel'),    cfg.artfctdef.tms.channel     = 'all';       end
    if ~isfield(cfg.artfctdef.tms,'trlpadding'), cfg.artfctdef.tms.trlpadding  = 0.1;         end
    if ~isfield(cfg.artfctdef.tms,'fltpadding'), cfg.artfctdef.tms.fltpadding  = 0.1;         end
    if ~isfield(cfg.artfctdef.tms,'artpadding'), cfg.artfctdef.tms.artpadding  = 0.01;        end
    if ~isfield(cfg.artfctdef.tms,'cutoff'),     cfg.artfctdef.tms.cutoff      = 4;           end
    % construct a temporary configuration that can be passed onto artifact_zvalue
    tmpcfg                  = [];
    tmpcfg.trl              = cfg.trl;
    tmpcfg.artfctdef.zvalue = cfg.artfctdef.tms;
    if isfield(cfg, 'continuous'),   tmpcfg.continuous       = cfg.continuous;    end
    if isfield(cfg, 'dataformat'),   tmpcfg.dataformat       = cfg.dataformat;    end
    if isfield(cfg, 'headerformat'), tmpcfg.headerformat     = cfg.headerformat;  end
    % call the zvalue artifact detection function
    
    % the data is either passed into the function by the user or read from file with cfg.inputfile
    hasdata = exist('data', 'var');
    
    if hasdata
      % read the header
      cfg = ft_checkconfig(cfg, 'forbidden', {'dataset', 'headerfile', 'datafile'});
      fsample = data.fsample;
      [tmpcfg, artifact] = ft_artifact_zvalue(tmpcfg, data);
    else
      cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
      cfg = ft_checkconfig(cfg, 'required', {'headerfile', 'datafile'});
      hdr = ft_read_header(cfg.headerfile);
      fsample = hdr.Fs;
      tmpcfg.datafile    = cfg.datafile;
      tmpcfg.headerfile  = cfg.headerfile;
      [tmpcfg, artifact] = ft_artifact_zvalue(tmpcfg);
    end
    cfg.artfctdef.tms = tmpcfg.artfctdef.zvalue;
    
    % adjust artifact definition so that Nx2 matrix contains detected TMS
    % events with user-specified pre- and post stimulus period included.
    % The reason for this is that ft_artifact_zvalue centers the period
    % marked as artifactual around the detected event. In the case of a TMS
    % pulse the window you would like to mark as artifactual is not
    % symmetrical around the onset of the pulse.
    
    % get values and express in samples
    prestim = round(cfg.prestim * fsample);
    poststim = round(cfg.poststim * fsample);
    
    % adjust Nx2 artifact matrix to be centered non-symmetrically around
    % detected TMS-pulse
    artifact(:,1) = (artifact(:,1)+artifact(:,2))./2 - prestim;
    artifact(:,2) = artifact(:,1) + poststim;
    cfg.artfctdef.tms.artifact = artifact;
    
  case 'marker'
    % Check if the cfg is correct for this method
    cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
    ft_checkconfig(cfg, 'required','trialdef');
    cfg.trialfun = ft_getopt(cfg, 'trialfun', 'ft_trialfun_general');
    trialdef = cfg.trialdef;
    trialdef.prestim = cfg.prestim;
    trialdef.poststim = cfg.poststim;
    cfg.trialdef = ft_checkconfig(trialdef,'required',{'eventvalue','eventtype'});
    
    % Get the trialfun
    cfg.trialfun = ft_getuserfun(cfg.trialfun, 'trialfun');
    
    % Evaluate the trialfun
    fprintf('evaluating trialfunction ''%s''\n', func2str(cfg.trialfun));
    trl   = feval(cfg.trialfun, cfg);
    
    % Prepare the found events for output
    artifact = trl(:,1:2);
    cfg.artfctdef.tms.artifact = artifact;
    fprintf('found %d events\n', size(artifact,1));
  otherwise
    error('unsupported method'); % This should be redundant as ft_checkconfig does not allow other methods than the supported ones.
end

cfg = rmfield(cfg, 'method'); % FIXME - not removing this causes problems when passing to ft_preprocessing

