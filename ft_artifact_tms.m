function [cfg, artifact] = ft_artifact_tms(cfg, data)

% FT_ARTIFACT_TMS reads the data segments of interest from file and identifies
% artefacts in EEG recordings that were done during TMS stimulation.
%
% Use as
%   [cfg, artifact] = ft_artifact_tms(cfg)
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
%   [cfg, artifact] = ft_artifact_tms(cfg, data)
% where the input data is a structure as obtained from FT_PREPROCESSING.
%
% In both cases the configuration should also contain
%   cfg.trl         = structure that defines the data segments of interest, see FT_DEFINETRIAL
%   cfg.continuous  = 'yes' or 'no' whether the file contains continuous data (default = 'yes')
% and
%   cfg.method      = 'detect' or 'marker', see below.
%   cfg.prestim     = scalar, time in seconds prior to onset of detected event to mark as artifactual (default = 0.005 seconds)
%   cfg.poststim    = scalar, time in seconds post onset of detected even to mark as artifactual (default = 0.010 seconds)
%
% The different methods are described in detail below.
%
% With cfg.method='detect', TMS-artifact are detected on basis of transient
% high-amplidude gradients that are typical for TMS-pulses. The data is preprocessed
% (again) with the following settings, which are optimal for identifying TMS-pulses.
% Artifacts are identified by means of thresholding the z-transformed value of the
% preprocessed data. This method acts as a wrapper around FT_ARTIFACT_ZVALUE.
%   cfg.artfctdef.tms.derivative  = 'yes'
%   cfg.artfctdef.tms.channel     = Nx1 cell-array with selection of channels, see FT_CHANNELSELECTION for details
%   cfg.artfctdef.tms.cutoff      = z-value at which to threshold (default = 4)
%   cfg.artfctdef.tms.trlpadding  = 0.1
%   cfg.artfctdef.tms.fltpadding  = 0.1
%   cfg.artfctdef.tms.artpadding  = 0.01
% Be aware that if one artifact falls within this specified range of another
% artifact, both artifact will be counted as one. Depending on cfg.prestim and
% cfg.poststim you may not mark enough data as artifactual.
%
% With cfg.method='marker', TMS-artifact onsets and offsets are based on
% markers/triggers that are written into the EEG dataset. This method acts as a
% wrapper around FT_DEFINETRIAL to determine on- and offsets of TMS pulses by reading
% markers in the EEG.
%   cfg.trialfun            = function name, see below (default = 'ft_trialfun_general')
%   cfg.trialdef.eventtype  = 'string'
%   cfg.trialdef.eventvalue = number, string or list with numbers or strings
% The cfg.trialfun option is a string containing the name of a function that you
% wrote yourself and that FT_ARTIFACT_TMS will call. The function should take the
% cfg-structure as input and should give a NxM matrix with M>=3 in the same format as
% "trl" as the output. You can add extra custom fields to the configuration structure
% to pass as arguments to your own trialfun. Furthermore, inside the trialfun you can
% use the FT_READ_EVENT function to get the event information from your data file.
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

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});
cfg = ft_checkconfig(cfg, 'required', 'method');
cfg = ft_checkconfig(cfg, 'allowedval', {'method', 'detect', 'marker'});

% set default rejection parameters
if ~isfield(cfg, 'artfctdef'),                       cfg.artfctdef                   = [];        end
if ~isfield(cfg.artfctdef, 'tms'),                   cfg.artfctdef.tms               = [];        end
if ~isfield(cfg, 'method'),                          cfg.method                      = 'detect';  end
if ~isfield(cfg, 'prestim'),                         cfg.prestim                     = 0.005;     end
if ~isfield(cfg, 'poststim'),                        cfg.poststim                    = 0.010;     end

if isfield(cfg.artfctdef.tms, 'artifact')
  fprintf('tms artifact detection has already been done, retaining artifacts\n');
  artifact = cfg.artfctdef.tms.artifact;
  return
end

switch cfg.method
  case 'detect'
    % settings for preprocessing
    if ~isfield(cfg.artfctdef.tms, 'derivative'),   cfg.artfctdef.tms.derivative   = 'yes';    end
    % settings for the zvalue subfunction
    if ~isfield(cfg.artfctdef.tms, 'method'),     cfg.artfctdef.tms.method      = 'zvalue';    end
    if ~isfield(cfg.artfctdef.tms, 'channel'),    cfg.artfctdef.tms.channel     = 'all';       end
    if ~isfield(cfg.artfctdef.tms, 'trlpadding'), cfg.artfctdef.tms.trlpadding  = 0.1;         end
    if ~isfield(cfg.artfctdef.tms, 'fltpadding'), cfg.artfctdef.tms.fltpadding  = 0.1;         end
    if ~isfield(cfg.artfctdef.tms, 'artpadding'), cfg.artfctdef.tms.artpadding  = 0.01;        end
    if ~isfield(cfg.artfctdef.tms, 'cutoff'),     cfg.artfctdef.tms.cutoff      = 4;           end
    
    % the data is either passed into the function by the user or read from file with cfg.inputfile
    hasdata = exist('data', 'var');
    
    % construct a temporary configuration that can be passed onto artifact_zvalue
    tmpcfg = keepfields(cfg, {'trl', 'continuous', 'dataset', 'datafile', 'headerfile', 'dataformat', 'headerformat'});
    tmpcfg.artfctdef.zvalue = cfg.artfctdef.tms;
    
    if hasdata
      fsample = data.fsample; % get the sampling rate from the data structure
      [tmpcfg, artifact] = ft_artifact_zvalue(tmpcfg, data);
    else
      [tmpcfg, artifact] = ft_artifact_zvalue(tmpcfg);
      % FT_ARTIFACT_ZVALUE will call FT_CHECKCONFIG with the dataset2files option
      hdr = ft_read_header(tmpcfg.headerfile, 'headerformat', tmpcfg.headerformat);
      fsample = hdr.Fs; % get the sampling rate from the file on disk
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
    ft_checkconfig(cfg, 'required', 'trialdef');
    cfg.trialfun = ft_getopt(cfg, 'trialfun', 'ft_trialfun_general');
    trialdef = cfg.trialdef;
    trialdef.prestim = cfg.prestim;
    trialdef.poststim = cfg.poststim;
    cfg.trialdef = ft_checkconfig(trialdef, 'required', {'eventvalue', 'eventtype'});
    
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
    % This should be redundant as ft_checkconfig does not allow other methods than the supported ones.
    ft_error('unsupported method');
end

% FIXME - not removing this causes problems when passing to ft_preprocessing
cfg = rmfield(cfg, 'method');
