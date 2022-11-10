function [cfg] = ft_databrowser(cfg, data)

% FT_DATABROWSER can be used for visual inspection of data. Artifacts that were
% detected by artifact functions (see FT_ARTIFACT_xxx functions where xxx is the type
% of artifact) are marked. Additionally data pieces can be marked and unmarked as
% artifact by manual selection. The output cfg contains the updated specification of
% the artifacts.
%
% Use as
%   [cfg] = ft_databrowser(cfg)
%   [cfg] = ft_databrowser(cfg, data)
% If you only specify the configuration structure, it should contain the name of the
% dataset on your hard disk (see below). If you specify input data, it should be a
% data structure as obtained from FT_PREPROCESSING or from FT_COMPONENTANALYSIS.
%
% If you want to browse data that is on disk, you have to specify
%   cfg.dataset                 = string with the filename
% Instead of specifying the dataset, you can also explicitely specify the name of the
% file containing the header information and the name of the file containing the
% data, using
%   cfg.datafile                = string with the filename
%   cfg.headerfile              = string with the filename
%
% The following configuration options are supported:
%   cfg.ylim                    = vertical scaling, can be 'maxmin', 'maxabs' or [ymin ymax] (default = 'maxabs')
%   cfg.zlim                    = color scaling to apply to component topographies, 'minmax', 'maxabs' (default = 'maxmin')
%   cfg.blocksize               = duration in seconds for cutting continuous data in segments
%   cfg.trl                     = structure that defines the data segments of interest, only applicable for trial-based data
%   cfg.continuous              = 'yes' or 'no', whether the data should be interpreted as continuous or trial-based
%   cfg.allowoverlap            = 'yes' or 'no', whether data that is overlapping in multiple trials is allowed (default = 'no')
%   cfg.channel                 = cell-array with channel labels, see FT_CHANNELSELECTION
%   cfg.channelclamped          = cell-array with channel labels, that when using the 'vertical' viewmode will always be shown at the bottom. This is useful for showing ECG/EOG channels along with the other channels
%   cfg.compscale               = string, 'global' or 'local', defines whether the colormap for the topographic scaling is applied per topography or on all visualized components (default = 'local')
%   cfg.viewmode                = string, 'vertical', 'butterfly', or 'component' for visualizing ICA/PCA topographies together with the timecourses (default = 'vertical')
%   cfg.plotlabels              = 'yes', 'no' or 'some', whether to plot channel labels in vertical viewmode. The option 'some' plots one label for every ten channels, which is useful if there are many channels (default = 'some')
%   cfg.plotevents              = 'no' or 'yes', whether to plot event markers (default = 'yes')
%   cfg.ploteventlabels         = 'type=value', 'colorvalue' (default = 'type=value')
%   cfg.artfctdef.xxx.artifact  = Nx2 matrix with artifact segments see FT_ARTIFACT_xxx functions
%   cfg.selectfeature           = string, name of feature to be selected/added (default = 'visual')
%   cfg.selectmode              = 'markartifact', 'markpeakevent', 'marktroughevent' (default = 'markartifact')
%   cfg.colorgroups             = 'sequential', 'allblack', 'labelcharN' (N = Nth character in label), 'chantype' or a vector with the length of the number of channels defining the groups (default = 'sequential')
%   cfg.linecolor               = string with line colors or Nx3 color map (default = customized lines map with 15 colors)
%   cfg.linewidth               = linewidth in points (default = 0.5)
%   cfg.linestyle               = linestyle/marker type, see options of the PLOT function (default = '-')
%   cfg.verticalpadding         = number or 'auto', padding to be added to top and bottom of plot to avoid channels largely dissappearing when viewmode = 'vertical'/'component'  (default = 'auto'). The padding is expressed as a proportion of the total height added to the top and bottom. The setting 'auto' determines the padding depending on the number of channels that are being plotted.
%   cfg.selfun                  = string, name of function that is evaluated using the right-click context menu. The selected data and cfg.selcfg are passed on to this function.
%   cfg.selcfg                  = configuration options for function in cfg.selfun
%   cfg.seldat                  = 'selected' or 'all', specifies whether only the currently selected or all channels will be passed to the selfun (default = 'selected')
%   cfg.visible                 = string, 'on' or 'off' whether figure will be visible (default = 'on')
%   cfg.position                = location and size of the figure, specified as a vector of the form [left bottom width height]
%   cfg.renderer                = string, 'opengl', 'zbuffer', 'painters', see MATLAB Figure Properties. If this function crashes, you should try 'painters'.
%   cfg.colormap                = string, or Nx3 matrix, see FT_COLORMAP
%
% The following options for the scaling of the EEG, EOG, ECG, EMG, MEG and NIRS channels
% is optional and can be used to bring the absolute numbers of the different
% channel types in the same range (e.g. fT and uV). The channel types are determined
% from the input data using FT_CHANNELSELECTION.
%   cfg.eegscale                = number, scaling to apply to the EEG channels prior to display
%   cfg.eogscale                = number, scaling to apply to the EOG channels prior to display
%   cfg.ecgscale                = number, scaling to apply to the ECG channels prior to display
%   cfg.emgscale                = number, scaling to apply to the EMG channels prior to display
%   cfg.megscale                = number, scaling to apply to the MEG channels prior to display
%   cfg.gradscale               = number, scaling to apply to the MEG gradiometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.magscale                = number, scaling to apply to the MEG magnetometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.nirsscale               = number, scaling to apply to the NIRS channels prior to display
%   cfg.mychanscale             = number, scaling to apply to the channels specified in cfg.mychan
%   cfg.mychan                  = Nx1 cell-array with selection of channels
%   cfg.chanscale               = Nx1 vector with scaling factors, one per channel specified in cfg.channel
%
% You can specify preprocessing options that are to be applied to the  data prior to
% display. Most options from FT_PREPROCESSING are supported. They should be specified
% in the sub-structure cfg.preproc like these examples
%   cfg.preproc.lpfilter        = 'no' or 'yes'  lowpass filter (default = 'no')
%   cfg.preproc.lpfreq          = lowpass  frequency in Hz
%   cfg.preproc.demean          = 'no' or 'yes', whether to apply baseline correction (default = 'no')
%   cfg.preproc.detrend         = 'no' or 'yes', remove linear trend from the data (done per trial) (default = 'no')
%   cfg.preproc.baselinewindow  = [begin end] in seconds, the default is the complete trial (default = 'all')
%
% In case of component viewmode, a layout is required. If no layout is specified, an
% attempt is made to construct one from the sensor definition that is present in the
% data or specified in the configuration.
%   cfg.layout                  = filename of the layout, see FT_PREPARE_LAYOUT
%   cfg.elec                    = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.grad                    = structure with gradiometer definition or filename, see FT_READ_SENS
% Additional plotting options for the component viewmode:
%   cfg.gridscale               = scalar, number of points along both directions for interpolation (default = 45 here)
%   cfg.shading                 = string, 'none', 'flat', 'interp' (default = 'flat')
%   cfg.interplimits            = string, 'electrodes' or 'mask' (default here = 'mask')
%   cfg.interpolation           = string, 'nearest', 'linear', 'natural', 'cubic' or 'v4' (default = 'v4')
%   cfg.contournum              = topoplot contour lines
%
% The default font size might be too small or too large, depending on the number of
% channels. You can use the following options to change the size of text inside the
% figure and along the axes.
%   cfg.fontsize                = number, fontsize inside the figure (default = 0.03)
%   cfg.fontunits               = string, can be 'normalized', 'points', 'pixels', 'inches' or 'centimeters' (default = 'normalized')
%   cfg.axisfontsize            = number, fontsize along the axes (default = 10)
%   cfg.axisfontunits           = string, can be 'normalized', 'points', 'pixels', 'inches' or 'centimeters' (default = 'points')
%
% When visually selection data, a right-click will bring up a context-menu containing
% functions to be executed on the selected data. You can use your own function using
% cfg.selfun and cfg.selcfg. You can use multiple functions by giving the names/cfgs
% as a cell-array.
%
% In butterfly and vertical mode, you can use the "identify" button to reveal the name of a
% channel. Please be aware that it searches only vertically. This means that it will
% return the channel with the amplitude closest to the point you have clicked at the
% specific time point. This might be counterintuitive at first.
%
% The "cfg.artfctdef" structure in the output cfg is comparable to the configuration
% used by the artifact detection functions like FT_ARTIFACT_ZVALUE and in
% FT_REJECTARTIFACT. It contains for each artifact type an Nx2 matrix in which the
% first column corresponds to the begin samples of an artifact period, the second
% column contains the end samples of the artifact periods.
%
% In case the databrowser crashes and you cannot close the window, use delete(gcf) to
% get rid of the figure.
%
% See also FT_PREPROCESSING, FT_REJECTARTIFACT, FT_ARTIFACT_EOG, FT_ARTIFACT_MUSCLE,
% FT_ARTIFACT_JUMP, FT_ARTIFACT_MANUAL, FT_ARTIFACT_THRESHOLD, FT_ARTIFACT_CLIP,
% FT_ARTIFACT_ECG, FT_COMPONENTANALYSIS

% Copyright (C) 2009-2022, Robert Oostenveld, Ingrid Nieuwenhuis
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
ft_preamble loadvar data
ft_preamble provenance data

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% the data can be passed as input arguments or can be read from disk
hasdata = exist('data', 'var');
hascomp = hasdata && ft_datatype(data, 'comp'); % can be 'raw+comp' or 'timelock+comp'

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'unused',     {'comps', 'inputfile', 'outputfile'});
cfg = ft_checkconfig(cfg, 'renamed',    {'zscale', 'ylim'});
cfg = ft_checkconfig(cfg, 'renamedval', {'ylim', 'auto', 'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval', {'selectmode', 'mark', 'markartifact'});
cfg = ft_checkconfig(cfg, 'renamedval', {'ploteventlabels', 'colorvalue', 'value'});
cfg = ft_checkconfig(cfg, 'renamed',    {'elecfile', 'elec'});
cfg = ft_checkconfig(cfg, 'renamed',    {'gradfile', 'grad'});
cfg = ft_checkconfig(cfg, 'renamed',    {'optofile', 'opto'});
cfg = ft_checkconfig(cfg, 'renamed',    {'channelcolormap', 'linecolor'});
cfg = ft_checkconfig(cfg, 'renamed',    {'anonimize', 'anonymize'}); % fix typo in previous version of the code
cfg = ft_checkconfig(cfg, 'renamed',    {'anonymise', 'anonymize'}); % use North American and Oxford British spelling
cfg = ft_checkconfig(cfg, 'renamed',    {'newfigure', 'figure'});
cfg = ft_checkconfig(cfg, 'deprecated', {'selectfeature'}); % please specify cfg.artfctdef.xxx and cfg.artfctdef.yyy for each feature

% ensure that the preproc specific options are located in the cfg.preproc substructure
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});

% set the defaults
cfg.ylim                = ft_getopt(cfg, 'ylim', 'maxabs');
cfg.zlim                = ft_getopt(cfg, 'zlim', 'maxmin');
cfg.artfctdef           = ft_getopt(cfg, 'artfctdef', struct);
cfg.selectmode          = ft_getopt(cfg, 'selectmode', 'markartifact');
cfg.blocksize           = ft_getopt(cfg, 'blocksize');                 % now used for both continuous and non-continuous data, defaulting done below
cfg.preproc             = ft_getopt(cfg, 'preproc');                   % see preproc for options
cfg.selfun              = ft_getopt(cfg, 'selfun');                    % default functions are 'simpleFFT', 'multiplotER', 'topoplotER', 'topoplotVAR', 'movieplotER'
cfg.selcfg              = ft_getopt(cfg, 'selcfg');                    % defaulting done below, requires layouts/etc to be processed
cfg.seldat              = ft_getopt(cfg, 'seldat', 'current');
cfg.colorgroups         = ft_getopt(cfg, 'colorgroups', 'sequential');
cfg.linecolor           = ft_getopt(cfg, 'linecolor', []); % the default is defined in lineattributes_common
cfg.linestyle           = ft_getopt(cfg, 'linestyle', '-');
cfg.linewidth           = ft_getopt(cfg, 'linewidth', 0.5);
cfg.eegscale            = ft_getopt(cfg, 'eegscale');
cfg.eogscale            = ft_getopt(cfg, 'eogscale');
cfg.ecgscale            = ft_getopt(cfg, 'ecgscale');
cfg.emgscale            = ft_getopt(cfg, 'emgscale');
cfg.megscale            = ft_getopt(cfg, 'megscale');
cfg.magscale            = ft_getopt(cfg, 'magscale');
cfg.gradscale           = ft_getopt(cfg, 'gradscale');
cfg.chanscale           = ft_getopt(cfg, 'chanscale');
cfg.mychanscale         = ft_getopt(cfg, 'mychanscale');
cfg.mychan              = ft_getopt(cfg, 'mychan');
cfg.layout              = ft_getopt(cfg, 'layout');
cfg.colormap            = ft_getopt(cfg, 'colormap', 'default');
cfg.plotlabels          = ft_getopt(cfg, 'plotlabels', 'some');
cfg.event               = ft_getopt(cfg, 'event');                       % this only exists for backward compatibility and should not be documented
cfg.continuous          = ft_getopt(cfg, 'continuous');                  % the default is set further down in the code, conditional on the input data
cfg.precision           = ft_getopt(cfg, 'precision', 'double');
cfg.compscale           = ft_getopt(cfg, 'compscale', 'local');
cfg.renderer            = ft_getopt(cfg, 'renderer');
cfg.fontsize            = ft_getopt(cfg, 'fontsize', 12);
cfg.fontunits           = ft_getopt(cfg, 'fontunits', 'points');         % inches, centimeters, normalized, points, pixels
cfg.editfontsize        = ft_getopt(cfg, 'editfontsize', 12);
cfg.editfontunits       = ft_getopt(cfg, 'editfontunits', 'points');     % inches, centimeters, normalized, points, pixels
cfg.axisfontsize        = ft_getopt(cfg, 'axisfontsize', 10);
cfg.axisfontunits       = ft_getopt(cfg, 'axisfontunits', 'points');     % inches, centimeters, normalized, points, pixels
cfg.verticalpadding     = ft_getopt(cfg, 'verticalpadding', 'auto');
cfg.allowoverlap        = ft_getopt(cfg, 'allowoverlap', 'no');          % for ft_fetch_data
cfg.contournum          = ft_getopt(cfg, 'contournum', 0);               % topoplot contour lines
cfg.trl                 = ft_getopt(cfg, 'trl');
cfg.gridscale           = ft_getopt(cfg, 'gridscale', 45);
cfg.shading             = ft_getopt(cfg, 'shading', 'flat');
cfg.interplimits        = ft_getopt(cfg, 'interplim', 'mask');
cfg.interpolation       = ft_getopt(cfg, 'interpmethod', 'v4');
cfg.channelclamped      = ft_getopt(cfg, 'channelclamped');
% set the defaults for plotting the events
cfg.plotevents          = ft_getopt(cfg, 'plotevents', 'yes');
cfg.ploteventlabels     = ft_getopt(cfg, 'ploteventlabels', 'type=value');
cfg.eventalpha          = ft_getopt(cfg, 'artifactalpha', 0.2);          % for the opacity of events
% set the defaults for plotting the artifacts
cfg.plotartifacts       = ft_getopt(cfg, 'plotartifacts', 'yes');
cfg.plotartifactlabels  = ft_getopt(cfg, 'plotartifactlabels', '');
cfg.artifactalpha       = ft_getopt(cfg, 'artifactalpha', 0.2);          % for the opacity of artifacts
% add some defaults for preprocessing, none of these is active but they will show up with the cfg.preproc button in the user interface
cfg.preproc.demean      = ft_getopt(cfg.preproc, 'demean', 'no');
cfg.preproc.lpfilter    = ft_getopt(cfg.preproc, 'lpfilter', 'no');
cfg.preproc.lpfreq      = ft_getopt(cfg.preproc, 'lpfreq', 30);
cfg.preproc.hpfilter    = ft_getopt(cfg.preproc, 'hpfilter', 'no');
cfg.preproc.hpfreq      = ft_getopt(cfg.preproc, 'hpfreq', 0.5);

% construct the low-level options as key-value pairs, these are passed to FT_READ_HEADER
headeropt = {};
headeropt  = ft_setopt(headeropt, 'headerformat',   ft_getopt(cfg, 'headerformat'));        % is passed to low-level function, empty implies autodetection
headeropt  = ft_setopt(headeropt, 'readbids',       ft_getopt(cfg, 'readbids'));            % is passed to low-level function
headeropt  = ft_setopt(headeropt, 'coordsys',       ft_getopt(cfg, 'coordsys', 'head'));    % is passed to low-level function
headeropt  = ft_setopt(headeropt, 'coilaccuracy',   ft_getopt(cfg, 'coilaccuracy'));        % is passed to low-level function
headeropt  = ft_setopt(headeropt, 'checkmaxfilter', ft_getopt(cfg, 'checkmaxfilter'));      % this allows to read non-maxfiltered neuromag data recorded with internal active shielding
headeropt  = ft_setopt(headeropt, 'chantype',       ft_getopt(cfg, 'chantype', {}));        % 2017.10.10 AB required for NeuroOmega files

% construct the low-level options as key-value pairs, these are passed to FT_READ_EVENT
eventopt = {};
eventopt  = ft_setopt(eventopt, 'eventformat',   ft_getopt(cfg, 'eventformat'));        % is passed to low-level function, empty implies autodetection
eventopt  = ft_setopt(eventopt, 'readbids',      ft_getopt(cfg, 'readbids'));           % is passed to low-level function

if isempty(ft_getopt(cfg, 'viewmode'))
  % can be 'butterfly', 'vertical', or 'component'
  if hascomp
    cfg.viewmode = 'component';
  else
    cfg.viewmode = 'vertical';
  end
end

if isempty(ft_getopt(cfg, 'colorgroups'))
  % can be 'sequential', 'allblack', 'labelcharN', 'chantype', or a vector with the length of the number of channels defining the groups
  if hascomp
    cfg.colorgroups = 'allblack';
  else
    cfg.colorgroups = 'sequential';
  end
end

if ~isempty(cfg.chanscale)
  if ~isfield(cfg, 'channel')
    ft_warning('ignoring cfg.chanscale; this should only be used when an explicit channel selection is being made');
    cfg.chanscale = [];
  elseif numel(cfg.channel) ~= numel(cfg.chanscale)
    ft_error('cfg.chanscale should have the same number of elements as cfg.channel');
  end
  % make sure chanscale is a column vector, not a row vector
  cfg.chanscale = cfg.chanscale(:);
end

if ~isempty(cfg.mychanscale) && ~isfield(cfg, 'mychan')
  ft_warning('ignoring cfg.mychanscale; no channels specified in cfg.mychan');
  cfg.mychanscale = [];
end

if isempty(ft_getopt(cfg, 'channel'))
  if hascomp
    if size(data.topo,2)>9
      cfg.channel = 1:10;
    else
      cfg.channel = 1:size(data.topo,2);
    end
  else
    cfg.channel = 'all';
  end
end

if strcmp(cfg.viewmode, 'component')
  % read or create the topographic layout that will be used for the topoplots
  tmpcfg = keepfields(cfg, {'layout', 'rows', 'columns', 'commentpos', 'skipcomnt', 'scalepos', 'skipscale', 'projection', 'viewpoint', 'rotate', 'width', 'height', 'elec', 'grad', 'opto', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
  if hasdata
    % select those channels from the layout that are relevant for the
    % decomposition that is being plotted
    tmpcfg.channel = data.topolabel;
    cfg.layout = ft_prepare_layout(tmpcfg, data);
  else
    cfg.layout = ft_prepare_layout(tmpcfg);
  end
end

if isempty(fieldnames(cfg.artfctdef)) % note that isempty(struct()) returns false
  % by default allow the user to specify visual artifacts
  cfg.artfctdef.visual.artifact = zeros(0,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the defaults and do some preparation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if hasdata
  % save whether data came from a timelock structure
  istimelock = strcmp(ft_datatype(data), 'timelock');

  if ~isempty(cfg.trl)
    % reduce the input data to the requested trial segments
    tmpcfg = keepfields(cfg, {'trl'});
    data = ft_redefinetrial(tmpcfg, data);
    [cfg, data] = rollback_provenance(cfg, data);
  end

  % check if the input data is valid for this function
  data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');
  % fetch the header from the data structure in memory
  hdr = ft_fetch_header(data);

  if isfield(data, 'cfg') && ~isempty(ft_findcfg(data.cfg, 'origfs'))
    % don't use the events in case the data has been resampled
    ft_warning('the data has been resampled, not showing the events');
    event = [];
  elseif istimelock
    % don't use the events in case the data has been averaged
    ft_warning('the data has been averaged, not showing the events');
    event = [];
  elseif ~isempty(cfg.event)
    % use the events that the user passed in the configuration
    event = cfg.event;
  else
    % fetch the events from the data structure in memory
    event = ft_fetch_event(data);
  end

  cfg.channel = ft_channelselection(cfg.channel, hdr.label);
  chansel = match_str(data.label, cfg.channel);
  Nchans  = length(chansel);

  if isempty(cfg.continuous)
    if numel(data.trial) == 1 && ~istimelock
      cfg.continuous = 'yes';
    else
      cfg.continuous = 'no';
    end
  else
    if strcmp(cfg.continuous, 'yes') && (numel(data.trial) > 1)
      ft_warning('interpreting trial-based data as continuous, the time axis now corresponds to the continuous data with the first sample being t=0')
    end
  end

  % this is how the input data is segmented
  trlorg = sampleinfo2trl(data);

else
  % check if the input cfg is valid for this function
  cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
  cfg = ft_checkconfig(cfg, 'required', {'headerfile', 'datafile'});
  cfg = ft_checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
  cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});
  % read the header from file
  hdr = ft_read_header(cfg.headerfile, headeropt{:});

  if isempty(cfg.continuous)
    if hdr.nTrials==1
      cfg.continuous = 'yes';
    else
      cfg.continuous = 'no';
    end
  end

  if ~isempty(cfg.event)
    % use the events that the user passed in the configuration
    event = cfg.event;
  else
    % read the events from file
    event = ft_read_event(cfg.dataset, eventopt{:});
  end

  cfg.channel = ft_channelselection(cfg.channel, hdr.label);
  chansel = match_str(hdr.label, cfg.channel);
  Nchans  = length(chansel);

  if ~isfield(cfg, 'trl') || isempty(cfg.trl)
    % treat the data as continuous if possible, otherwise define all trials as indicated in the header
    if strcmp(cfg.continuous, 'yes')
      trlorg = zeros(1, 3);
      trlorg(1,1) = 1;
      trlorg(1,2) = hdr.nSamples*hdr.nTrials;
      trlorg(1,3) = -hdr.nSamplesPre;
    else
      trlorg = zeros(hdr.nTrials, 3);
      for i=1:hdr.nTrials
        trlorg(i,1) = (i-1)*hdr.nSamples + 1;
        trlorg(i,2) = (i  )*hdr.nSamples    ;
        trlorg(i,3) = -hdr.nSamplesPre;
      end
    end
  elseif ischar(cfg.trl)
    % load the trial information from file
    trlorg = loadvar(cfg.trl, 'trl');
  else
    trlorg = cfg.trl;
  end

end % if hasdata

% the code below expects an Nx3 matrix with begsample, endsample and offset
if istable(trlorg)
  trlorg = table2array(trlorg(:,1:3));
else
  trlorg = trlorg(:,1:3);
end

Ntrials = size(trlorg, 1);

if strcmp(cfg.continuous, 'no') && isempty(cfg.blocksize)
  cfg.blocksize = (trlorg(1,2) - trlorg(1,1)+1) ./ hdr.Fs;
elseif strcmp(cfg.continuous, 'yes') && isempty(cfg.blocksize)
  cfg.blocksize = 1;
end

if cfg.blocksize<round(10*1/hdr.Fs)
  ft_warning('the blocksize is very small given the samping rate, increasing blocksize to 10 samples');
  cfg.blocksize = round(10*1/hdr.Fs);
end

% FIXME make a check for the consistency of cfg.continuous, cfg.blocksize, cfg.trl and the data header

if Nchans == 0
  ft_error('no channels to display');
end

if Ntrials == 0
  ft_error('no trials to display');
end

% determine the vertical scaling
if ischar(cfg.ylim)
  if hasdata
    sel = 1;
    while all(isnan(reshape(data.trial{sel}(chansel,:),[],1)))
      sel = sel+1;
    end
    % the first trial is used to determine the vertical scaling
    dat = data.trial{sel}(chansel,:);
  else
    % read one second (or one block) of data to determine the vertical scaling
    begsample = 1;
    endsample = min(round(hdr.Fs), hdr.nSamples);
    dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chansel, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat, headeropt{:});
  end % if hasdata
  % convert the data to another numeric precision, i.e. double, single or int32
  if ~isempty(cfg.precision)
    dat = cast(dat, cfg.precision);
  end
  minval = min(dat(:));
  maxval = max(dat(:));
  switch cfg.ylim
    case 'maxabs'
      maxabs   = max(abs([minval maxval]));
      scalefac = 10^(fix(log10(maxabs)));
      if scalefac==0
        % this happens if the data is all zeros
        scalefac=1;
      end
      maxabs   = (round(maxabs / scalefac * 100) / 100) * scalefac;
      cfg.ylim = [-maxabs maxabs];
    case 'maxmin'
      if minval==maxval
        % this happens if the data is constant, e.g. all zero or clipping
        minval = minval - eps;
        maxval = maxval + eps;
      end
      cfg.ylim = [minval maxval];
    otherwise
      ft_error('unsupported value for cfg.ylim');
  end % switch ylim
  if strcmp(cfg.viewmode, 'vertical') || strcmp(cfg.viewmode, 'component')
    % it is OK to have some overlap between the traces in vertical and component viewmodes
    % but butterfly plots should stay within the min/max boundaries
    scale_adjust = 5;
    cfg.ylim = cfg.ylim/scale_adjust;
  end
else
  if (numel(cfg.ylim) ~= 2) || ~isnumeric(cfg.ylim)
    ft_error('cfg.ylim needs to be a 1x2 vector [ymin ymax], describing the upper and lower limits')
  end
end

% determine the coloring of channels
if hasdata
  linecolor = lineattributes_common(cfg, data);
else
  linecolor = lineattributes_common(cfg, hdr);
end

% collect the artifacts from cfg.artfctdef.xxx.artifact
artlabel = fieldnames(cfg.artfctdef);
sel      = zeros(size(artlabel));
artifact = cell(size(artlabel));
for i=1:length(artlabel)
  sel(i) = isfield(cfg.artfctdef.(artlabel{i}), 'artifact');
  if sel(i)
    artifact{i} = cfg.artfctdef.(artlabel{i}).artifact;
    ft_info('detected %3d %s artifacts\n', size(artifact{i}, 1), artlabel{i});
  end
end
% continue with the subset of artfctdef fields that actually contain artifacts
artifact = artifact(sel==1);
artlabel = artlabel(sel==1);

% make artdata representing all artifacts in a "raw data" format
endsample = max(trlorg(:,2));

artdata = [];
artdata.trial{1}       = artifact2boolvec(artifact, 'endsample', endsample); % every artifact is a "channel"
artdata.time{1}        = offset2time(0, hdr.Fs, endsample);
artdata.label          = artlabel;
artdata.fsample        = hdr.Fs;
artdata.cfg.trl        = [1 endsample 0];

% determine the unique artifact types and corresponding colors, this only needs to be done once
artifacttypes = artlabel;
artifactcolors = colorcheck([0.9686 0.7608 0.7686; 0.7529 0.7098 0.9647; 0.7373 0.9725 0.6824; 0.8118 0.8118 0.8118; 0.9725 0.6745 0.4784; 0.9765 0.9176 0.5686; 0.6863 1 1; 1 0.6863 1; 0 1 0.6000], numel(artdata.label));

% determine the unique event types and corresponding colors, this only needs to be done once
if ~isempty(event) && isstruct(event)
  eventtypes  = unique({event.type});
  eventcolors = colorcheck('krbgmcy', numel(eventtypes));
  % durations and offsets can be either empty or should be numeric values, see FT_READ_EVENT
  % the code further down expects them to be numeric values, so change them to zero
  for i=1:numel(event)
    if isempty(event(i).duration)
      event(i).duration = 0;
    end
    if isempty(event(i).offset)
      event(i).offset = 0;
    end
  end
else
  eventtypes = {};
  eventcolors = '';
end

ft_info('the different artifact types correspond to the following colors:');
tmp = pad(artifacttypes); % pad the artifacttypes to the same width
for i=1:length(artifacttypes)
  ft_info('  %s = %s\n', tmp{i}, htmlcolors(artifactcolors(i,:)));
end
clear tmp

ft_info('the different event types correspond to the following colors:');
tmp = pad(eventtypes); % pad the eventtypes to the same width
for i=1:length(eventtypes)
  ft_info('  %s = %s\n', tmp{i}, htmlcolors(eventcolors(i,:)));
end
clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up default functions to be available in the right-click menu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cfg.selfun - labels that are presented in rightclick menu, and is appended using ft_getuserfun(..., 'browse') later on to create a function handle
% cfg.selcfg - cfgs for functions to be executed

if ~isempty(cfg.selfun) || ~isempty(cfg.selcfg)
  if ischar(cfg.selfun)
    cfg.selfun = {cfg.selfun};
  end
  if isstruct(cfg.selcfg)
    cfg.selcfg = {cfg.selcfg};
  end
elseif isempty(cfg.selfun) && isempty(cfg.selcfg)
  % simplefft
  cfg.selcfg{1} = [];
  cfg.selcfg{1}.linecolor = linecolor;
  cfg.selfun{1} = 'simpleFFT';
  % multiplotER
  cfg.selcfg{2} = [];
  cfg.selcfg{2}.linecolor = linecolor;
  cfg.selcfg{2}.layout = cfg.layout;
  cfg.selcfg{2}.colorgroups = 'sequential';
  cfg.selfun{2} = 'multiplotER';
  % topoplotER
  cfg.selcfg{3} = [];
  cfg.selcfg{3}.linecolor = linecolor;
  cfg.selcfg{3}.layout    = cfg.layout;
  cfg.selcfg{3}.colormap  = cfg.colormap;
  cfg.selfun{3} = 'topoplotER';
  % topoplotVAR
  cfg.selcfg{4} = [];
  cfg.selcfg{4}.layout   = cfg.layout;
  cfg.selcfg{4}.colormap = cfg.colormap;
  cfg.selfun{4} = 'topoplotVAR';
  % movieplotER
  cfg.selcfg{5} = [];
  cfg.selcfg{5}.layout = cfg.layout;
  cfg.selcfg{5}.interactive = 'yes';
  cfg.selfun{5} = 'movieplotER';
  % audiovideo
  cfg.selcfg{6} = [];
  cfg.selcfg{6}.audiofile = ft_getopt(cfg, 'audiofile');
  cfg.selcfg{6}.videofile = ft_getopt(cfg, 'videofile');
  cfg.selcfg{6}.anonymize = ft_getopt(cfg, 'anonymize');
  cfg.selfun{6} = 'audiovideo';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the data structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the structure "opt" represents the global data+settings, it should contain
% - the original data, epoched or continuous
% - the artifacts represented as continuous data
% - the redraw_cb settings
% - the preproc   settings
% - the select_range_cb settings (also used in keypress_cb)

% these elements are stored inside the figure, so that the callback routines can modify them
opt = [];
if hasdata
  opt.orgdata   = data;
else
  opt.orgdata   = [];        % this means that it will read from cfg.dataset
  opt.headeropt = headeropt; % these are passed to FT_READ_HEADER and FT_READ_DATA
end
if strcmp(cfg.continuous, 'yes')
  opt.trialviewtype = 'segment';
else
  opt.trialviewtype = 'trial';
end
opt.hdr             = hdr;
opt.trlop           = 1;          % the active trial being displayed
opt.artsel          = 1;          % the currently selected artifact/feature
opt.trlorg          = trlorg;
opt.fsample         = hdr.Fs;
opt.linecolor       = linecolor;
opt.changedchanflg  = true;       % trigger for redrawing channel labels and preparing layout again (see bug 2065 and 2878)
opt.cleanup         = false;      % this is needed for a corrent handling if the figure is closed (either in the corner or by "q")
opt.chanindx        = [];         % this is used to check whether the component topographies need to be redrawn
opt.artdata         = artdata;
opt.artifacttypes   = artifacttypes;
opt.artifactcolors  = artifactcolors;
opt.event           = event;
opt.eventtypes      = eventtypes;
opt.eventcolors     = eventcolors;
opt.nanpaddata      = []; % this is used to allow horizontal scaling to be constant (when looking at last segment continuous data, or when looking at segmented/zoomed-out non-continuous data)
opt.trllock         = []; % this is used when zooming into trial based data

% save original layout when viewmode = component
if strcmp(cfg.viewmode, 'component')
  opt.layouttopo = cfg.layout;
end

% open a new figure with the specified settings
h = open_figure(keepfields(cfg, {'figure', 'position', 'visible', 'renderer'}));

% check if the colormap is in proper format and set it
if ~isempty(cfg.colormap)
  if ischar(cfg.colormap)
    cfg.colormap = ft_colormap(cfg.colormap);
  elseif iscell(cfg.colormap)
    cfg.colormap = ft_colormap(cfg.colormap{:});
  elseif isnumeric(cfg.colormap) && size(cfg.colormap,2)~=3
    ft_error('cfg.colormap must be Nx3');
  end
  set(h, 'colormap', cfg.colormap);
end

% put appdata in figure
setappdata(h, 'opt', opt);
setappdata(h, 'cfg', cfg);

% set interruptible to off, see bug 3123
set(h, 'Interruptible', 'off', 'BusyAction', 'queue'); % enforce busyaction to queue to be sure

% enable custom data cursor text
dcm = datacursormode(h);
set(dcm, 'updatefcn', @cb_datacursortext);

% set the figure window title
funcname = mfilename();
if ~hasdata
  if isfield(cfg, 'dataset')
    dataname = cfg.dataset;
  elseif isfield(cfg, 'datafile')
    dataname = cfg.datafile;
  else
    dataname = [];
  end
elseif isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)
  dataname = cfg.inputfile;
else
  dataname = inputname(2);
end
set(gcf, 'Name', sprintf('%d: %s: %s', double(gcf), funcname, join_str(', ',dataname)));
set(gcf, 'NumberTitle', 'off');

% set zoom option to on
% zoom(h, 'on')
% set(zoom(h), 'actionPostCallback', @zoom_drawlabels_cb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add the graphical user interface elements and callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these will be set using enable_callbacks after each redraw
set(h, 'KeyPressFcn',           []);
set(h, 'WindowButtonDownFcn',   []);
set(h, 'WindowButtonUpFcn',     []);
set(h, 'WindowButtonMotionFcn', []);

if ~isempty(findobj(h, 'tag', 'navigateui'))
  % assume that the graphical user interface elements are already present
  % this speeds up plotting in cases like this https://www.fieldtriptoolbox.org/example/video_eeg/

else
  % add the graphical user interface elements
  uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', opt.trialviewtype, 'userdata', 't')
  uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'leftarrow')
  uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'rightarrow')

  if strcmp(cfg.viewmode, 'component')
    uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'component', 'userdata', 'c')
  else
    uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'channel', 'userdata', 'c')
  end

  uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'uparrow')
  uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'downarrow')

  uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'horizontal', 'userdata', 'h')
  uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'shift+leftarrow')
  uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'shift+rightarrow')

  uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'vertical', 'userdata', 'v')
  uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'shift+downarrow')
  uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'shift+uparrow')

  % preproc button - allows updating the preprocessing options on the fly
  uicontrol('tag', 'rightcolumn', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'preproc', 'userdata', 'p', 'position', [0.91, 0.88, 0.08, 0.04])

  % identify button - to find label of nearest channel to datapoint
  uicontrol('tag', 'rightcolumn', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'identify', 'userdata', 'i', 'position', [0.91, 0.83, 0.08, 0.04])

  % viewmode button to toggle between vertical and butterfly
  if ~strcmp(cfg.viewmode, 'component')
    uicontrol('tag', 'rightcolumn', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'viewmode', 'userdata', 'm', 'position', [0.91, 0.78, 0.08, 0.04])
  end

  % legend artifacts/features
  for iArt = 1:length(artlabel)
    uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', artlabel{iArt}, 'userdata',  num2str(iArt),  'position', [0.91, 0.850 - ((2+iArt-1)*0.09), 0.08, 0.04], 'backgroundcolor', opt.artifactcolors(iArt,:))
    uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', ['shift+'   num2str(iArt)], 'position', [0.91, 0.805 - ((2+iArt-1)*0.09), 0.03, 0.04], 'backgroundcolor', opt.artifactcolors(iArt,:))
    uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', ['control+' num2str(iArt)], 'position', [0.96, 0.805 - ((2+iArt-1)*0.09), 0.03, 0.04], 'backgroundcolor', opt.artifactcolors(iArt,:))
  end
  if length(artlabel)>1 % highlight the first one as active
    arth = findobj(h, 'tag', 'artifactui');
    arth = arth(end:-1:1); % order is reversed so reverse it again
    hsel = [1 2 3] + (opt.artsel-1) .*3;
    set(arth(hsel), 'fontweight', 'bold')
  end

  ft_uilayout(h, 'tag', 'labels',  'width', 0.10, 'height', 0.05);
  ft_uilayout(h, 'tag', 'buttons', 'width', 0.05, 'height', 0.05);

  ft_uilayout(h, 'tag', 'labels',      'style', 'pushbutton', 'callback', @keypress_cb);
  ft_uilayout(h, 'tag', 'buttons',     'style', 'pushbutton', 'callback', @keypress_cb);
  ft_uilayout(h, 'tag', 'rightcolumn', 'style', 'pushbutton', 'callback', @keypress_cb);
  ft_uilayout(h, 'tag', 'artifactui',  'style', 'pushbutton', 'callback', @keypress_cb);

  ft_uilayout(h, 'tag', 'labels',  'retag', 'navigateui'); % this renames the existing tags
  ft_uilayout(h, 'tag', 'buttons', 'retag', 'navigateui'); % this renames the existing tags
  ft_uilayout(h, 'tag', 'navigateui', 'BackgroundColor', [0.8 0.8 0.8], 'hpos', 'auto', 'vpos', 0);

  % add a menu to the figure, but only if the current figure does not have subplots
  tmpcfg = cfg;
  if hasdata && isfield(data, 'cfg')
    tmpcfg.previous = data.cfg;
  end
  menu_fieldtrip(h, tmpcfg, false);

end % if findobj navigateui

definetrial_cb(h);
redraw_cb(h);
help_cb(h);

if nargout
  % wait until the user interface is closed, get the user data with the updated artifact details
  set(h, 'CloseRequestFcn', @cleanup_cb);

  while ishandle(h)
    uiwait(h);
    opt = getappdata(h, 'opt');
    if opt.cleanup
      delete(h);
    end
  end

  % add the updated artifact definitions to the output cfg
  for i=1:length(opt.artdata.label)
    cfg.artfctdef.(opt.artdata.label{i}).artifact = boolvec2artifact(opt.artdata.trial{1}(i,:));
  end

  % add the updated preproc to the output
  try
    browsecfg = getappdata(h, 'cfg');
    cfg.preproc = browsecfg.preproc;
  end

  % add the updated events to the output cfg
  cfg.event = opt.event;

end % if nargout

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous data
ft_postamble provenance

if ~nargout
  clear cfg
end

end % main function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end
end % function getparent


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cursortext = cb_datacursortext(h, eventdata)
pos = get(eventdata, 'Position');

linetype = getappdata(eventdata.Target, 'datacursor_linetype');

if strcmp(linetype, 'event')
  cursortext = sprintf('%s\nt = %g s', getappdata(eventdata.Target, 'datacursor_eventlabel'), getappdata(eventdata.Target, 'datacursor_eventtime'));

elseif strcmp(linetype, 'artifact')
  cursortext = sprintf('%s\nt = %g s', getappdata(eventdata.Target, 'datacursor_artifactlabel'), getappdata(eventdata.Target, 'datacursor_artifacttime'));

elseif strcmp(linetype, 'channel')
  % get plotted x axis
  plottedX = get(eventdata.Target, 'xdata');

  % determine values of data at real x axis
  timeAxis = getappdata(eventdata.Target, 'datacursor_xdata');
  dataAxis = getappdata(eventdata.Target, 'datacursor_ydata');
  tInd = nearest(plottedX, pos(1));

  % get label
  chanLabel = getappdata(eventdata.Target, 'datacursor_label');
  chanLabel = chanLabel{1};

  cursortext = sprintf('%s = %g\nt = %g', chanLabel, dataAxis(tInd), timeAxis(tInd));

else
  % explicitly tell the user there is no info because the x-axis and y-axis do not correspond to real data values (both are always between 0 and 1)
  cursortext = '<no cursor available>';
end
end % function cb_datacursortext


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION see also lineattributes_common
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function color = colorcheck(color, n)
% define the mapping between color characters and RGB values
name = 'rgbcmywk';
rgb = [
  1 0 0
  0 1 0
  0 0 1
  0 1 1
  1 0 1
  1 1 0
  1 1 1
  0 0 0
  ];
% ensure that all colors are represented as RGB
if ischar(color)
  original = color;
  color = nan(length(original), 3);
  for i=1:length(original)
    color(i,:) = rgb(name==original(i),:);
  end
end
% ensure that there are enough colors
color = repmat(color, ceil(n/size(color,1)), 1);
color = color(1:n,:);
end % function colorcheck


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanup_cb(h, eventdata)
opt = getappdata(h, 'opt');
opt.cleanup = true;
setappdata(h, 'opt', opt);
uiresume
end % function cleanup_cb


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function definetrial_cb(h, eventdata)
h = getparent(h);
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

if strcmp(cfg.continuous, 'no')

  % when zooming in, lock the trial! one can only go to the next trial when horizontal scaling doesn't segment the data - from ft-meeting: this might be relaxed later on - roevdmei
  if isempty(opt.trllock)
    opt.trllock = opt.trlop;
  end
  locktrllen = ((opt.trlorg(opt.trllock,2)-opt.trlorg(opt.trllock,1)+1) ./ opt.fsample);
  % if cfg.blocksize is close to the length of the locked trial, set it to that
  if (abs(locktrllen-cfg.blocksize) / locktrllen) < 0.1
    cfg.blocksize = locktrllen;
  end

  %%%%%%%%%
  % trial is locked, change subdivision of trial
  if cfg.blocksize < locktrllen
    % lock the trial if it wasn't locked (and thus trlop refers to the actual trial)
    if isempty(opt.trllock)
      opt.trllock = trlop;
    end
    % save current position if already
    if isfield(opt, 'trlvis')
      thissegbeg = opt.trlvis(opt.trlop,1);
    end
    datbegsample = min(opt.trlorg(opt.trllock,1));
    datendsample = max(opt.trlorg(opt.trllock,2));
    smpperseg  = round(opt.fsample * cfg.blocksize);
    begsamples = datbegsample:smpperseg:datendsample;
    endsamples = datbegsample+smpperseg-1:smpperseg:datendsample;
    offset     = (((1:numel(begsamples))-1)*smpperseg) + opt.trlorg(opt.trllock,3);
    if numel(endsamples)<numel(begsamples)
      endsamples(end+1) = datendsample;
    end
    trlvis = [];
    trlvis(:,1) = begsamples';
    trlvis(:,2) = endsamples';
    trlvis(:,3) = offset;
    % determine length of each trial, and determine the offset with the current requested zoom-level
    trllen   = (trlvis(:,2) - trlvis(:,1)+1);
    sizediff = smpperseg - trllen;
    opt.nanpaddata = sizediff;

    if isfield(opt, 'trlvis')
      % update the current trial counter and try to keep the current sample the same
      opt.trlop   = nearest(begsamples, thissegbeg);
    end
    % update trialviewtype
    opt.trialviewtype = 'trialsegment';
    % update button
    set(findobj(get(h, 'children'), 'string', 'trial'), 'string', 'segment');
    %%%%%%%%%


    %%%%%%%%%
    % trial is not locked, go to original trial division and zoom out
  elseif cfg.blocksize >= locktrllen
    trlvis = opt.trlorg;
    % set current trlop to locked trial if it was locked before
    if ~isempty(opt.trllock)
      opt.trlop = opt.trllock;
    end
    smpperseg  = round(opt.fsample * cfg.blocksize);
    % determine length of each trial, and determine the offset with the current requested zoom-level
    trllen   = (trlvis(:,2) - trlvis(:,1)+1);
    sizediff = smpperseg - trllen;
    opt.nanpaddata = sizediff;

    % update trialviewtype
    opt.trialviewtype = 'trial';
    % update button
    set(findobj(get(h, 'children'), 'string', 'trialsegment'), 'string', opt.trialviewtype);

    % release trial lock
    opt.trllock = [];
    %%%%%%%%%
  end

  % save trlvis
  opt.trlvis  = trlvis;

else
  % construct a trial definition for visualisation
  if isfield(opt, 'trlvis') % if present, remember where we were
    thistrlbeg = opt.trlvis(opt.trlop,1);
  end
  % look at cfg.blocksize and make opt.trl accordingly
  datbegsample = min(opt.trlorg(:,1));
  datendsample = max(opt.trlorg(:,2));
  smpperseg  = round(opt.fsample * cfg.blocksize);
  begsamples = datbegsample:smpperseg:datendsample;
  endsamples = datbegsample+smpperseg-1:smpperseg:datendsample;
  if numel(endsamples)<numel(begsamples)
    endsamples(end+1) = datendsample;
  end
  trlvis = [];
  trlvis(:,1) = begsamples';
  trlvis(:,2) = endsamples';
  % compute the offset. In case if opt.trlorg has multiple trials, the first sample is t=0, otherwise use the offset in opt.trlorg
  if size(opt.trlorg,1)==1
    offset = begsamples - repmat(begsamples(1), [1 numel(begsamples)]); % offset for all segments compared to the first
    offset = offset + opt.trlorg(1,3);
    trlvis(:,3) = offset;
  else
    offset = begsamples - repmat(begsamples(1), [1 numel(begsamples)]);
    trlvis(:,3) = offset;
  end

  if isfield(opt, 'trlvis')
    % update the current trial counter and try to keep the current sample the same
    % opt.trlop   = nearest(round((begsamples+endsamples)/2), thissample);
    opt.trlop   = nearest(begsamples, thistrlbeg);
  end
  opt.trlvis  = trlvis;

  % NaN-padding when horizontal scaling is bigger than the data
  % two possible situations, 1) zoomed out so far that all data is one segment, or 2) multiple segments but last segment is smaller than the rest
  sizediff = smpperseg-(endsamples-begsamples+1);
  opt.nanpaddata = sizediff;
end % if continuous
setappdata(h, 'opt', opt);
setappdata(h, 'cfg', cfg);
end % function definetrial_cb


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_range_cb(h, range, cmenulab)
h = getparent(h);
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

% the range should be in the displayed box
range(1) = max(opt.hpos-opt.width/2, range(1));
range(2) = max(opt.hpos-opt.width/2, range(2));
range(1) = min(opt.hpos+opt.width/2, range(1));
range(2) = min(opt.hpos+opt.width/2, range(2));
range = (range-(opt.hpos-opt.width/2)) / opt.width;          % left side of the box becomes 0, right side becomes 1
range = range * (opt.hlim(2) - opt.hlim(1)) + opt.hlim(1);   % 0 becomes hlim(1), 1 becomes hlim(2)

begsample = opt.trlvis(opt.trlop,1);
endsample = opt.trlvis(opt.trlop,2);
offset    = opt.trlvis(opt.trlop,3);

% determine the selection, it is now always based on begsample/endsample/offset -roevdmei
begsel = round(range(1)*opt.fsample+begsample-offset-1);
endsel = round(range(2)*opt.fsample+begsample-offset);

% the selection should always be confined to the current trial
begsel = max(begsample, begsel);
endsel = min(endsample, endsel);

% mark or execute selfun
if isempty(cmenulab)
  % the left button was clicked INSIDE a selected range, update the artifact definition or event

  if strcmp(cfg.selectmode, 'markartifact')
    % mark or unmark artifacts
    artval = opt.artdata.trial{1}(opt.artsel, begsel:endsel);
    artval = any(artval,1);
    if any(artval)
      fprintf('there is overlap with the active artifact (%s), disabling this artifact\n',opt.artdata.label{opt.artsel});
      opt.artdata.trial{1}(opt.artsel, begsel:endsel) = 0;
    else
      fprintf('there is no overlap with the active artifact (%s), marking this as a new artifact\n',opt.artdata.label{opt.artsel});
      opt.artdata.trial{1}(opt.artsel, begsel:endsel) = 1;
    end

    % redraw only when marking (so the focus doesn't go back to the databrowser after calling selfuns
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    redraw_cb(h);

  elseif strcmp(cfg.selectmode, 'markpeakevent') || strcmp(cfg.selectmode, 'marktroughevent')
    % mark or unmark events at the peak or trough in the selected window
    if ~isempty(opt.event) && any(intersect(begsel:endsel, [opt.event.sample]))
      fprintf('there is overlap with one or more event(s), removing this/these event(s)\n');
      ind_rem = intersect(begsel:endsel, [opt.event.sample]);
      for iRemove = 1:length(ind_rem)
        opt.event([opt.event.sample]==ind_rem(iRemove)) = [];
      end
    else
      dat = opt.curdata.trial{1}(:,(begsel-begsample+1):(endsel-begsample+1));
      if strcmp(cfg.selectmode, 'markpeakevent')
        fprintf('adding an event to the largest peak or maxium value\n');
        [val, indx] = max(dat(:));
        [chanindx, smpindx] = ind2sub(size(dat), indx);
        type = 'peak';
      elseif strcmp(cfg.selectmode, 'marktroughevent')
        fprintf('adding an event to the largest trough or negative value\n');
        [val, indx] = min(dat(:));
        [chanindx, smpindx] = ind2sub(size(dat), indx);
        type = 'trough';
      end
      event_new.type     = [type '_' opt.curdata.label{chanindx}];
      event_new.sample   = begsel + smpindx - 1;
      event_new.value    = val;
      event_new.duration = 1;
      event_new.offset   = 0;
      disp(event_new)
      if isempty(opt.event)
        % store the new event
        opt.event = event_new;
      else
        % append the new event
        opt.event(end+1) = event_new;
        % sort to keep the order consistent
        [dum, indx] = sort([opt.event.sample]);
        opt.event = opt.event(indx);
      end
    end % remove or add an event

    % redraw only when marking (so the focus doesn't go back to the databrowser after calling selfuns
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    redraw_cb(h);
  end

else
  % the right button was used to activate the context menu and the user made a selection from that menu
  % execute the corresponding function

  % get index into cfgs
  selfunind = strcmp(cfg.selfun, cmenulab);

  % cut out the requested data segment
  switch cfg.seldat
    case 'current'
      seldata             = keepfields(opt.curdata, {'label', 'grad', 'elec', 'opto', 'hdr'});
      seldata.trial{1}    = ft_fetch_data(opt.curdata, 'begsample', begsel, 'endsample', endsel);
    case 'all'
      seldata             = keepfields(opt.org, {'label', 'grad', 'elec', 'opto', 'hdr'});
      seldata.trial{1}    = ft_fetch_data(opt.orgdata, 'begsample', begsel, 'endsample', endsel);
  end
  seldata.time{1}     = offset2time(offset+begsel-begsample, opt.fsample, endsel-begsel+1);
  seldata.fsample     = opt.fsample;
  seldata.sampleinfo  = [begsel endsel];

  % prepare input
  funhandle = ft_getuserfun(cmenulab, 'browse');
  funcfg    = cfg.selcfg{selfunind};
  % construct a descriptive name for the new figure
  if ~strcmp(opt.trialviewtype, 'trialsegment')
    funcfg.figurename = sprintf('%s : %s %d/%d, time from %g to %g s', cmenulab, opt.trialviewtype, opt.trlop, size(opt.trlvis,1), seldata.time{1}(1), seldata.time{1}(end));
  else
    funcfg.figurename = sprintf('%s : trial %d/%d: segment: %d/%d , time from %g to %g s', cmenulab, opt.trllock, size(opt.trlorg,1), opt.trlop, size(opt.trlvis,1), seldata.time{1}(1), seldata.time{1}(end));
  end
  
  if ~isempty(opt.orgdata) && isfield(funcfg, 'linecolor')
    % the function that is executed does not know that only a subset of the channels will be passed in the input,
    % make sure that the funcfg's linecolor is consistent with this selection
    selchan = match_str(opt.orgdata.label, seldata.label);
    funcfg.linecolor = opt.linecolor(selchan, :);
  end

  feval(funhandle, funcfg, seldata);
end

end % function select_range_cb


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function preproc_cfg1_cb(h, eventdata)
parent = getparent(h);
cfg = getappdata(parent, 'cfg');

editfontsize = cfg.editfontsize;
editfontunits = cfg.editfontunits;

% parse cfg.preproc
if ~isempty(cfg.preproc)
  code = printstruct('cfg.preproc', cfg.preproc);
else
  code = '';
end

% add descriptive lines
code = {
  '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  '% Add or change options for on-the-fly preprocessing'
  '% Use as cfg.preproc.option=value'
  '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  code
  };

% make figure displaying the edit box
pph = figure;
axis off
% add save button
uicontrol('tag', 'preproccfg_l2', 'parent', pph, 'units', 'normalized', 'style', 'pushbutton', 'string', 'save and close', 'position', [0.81, 0.6 , 0.18, 0.10], 'callback', @preproc_cfg2_cb);

% add edit box
ppeh = uicontrol('style', 'edit');
set(pph, 'toolBar', 'none')
set(pph, 'menuBar', 'none')
set(pph, 'Name', 'cfg.preproc editor')
set(pph, 'NumberTitle', 'off')
set(ppeh, 'Units', 'normalized');
set(ppeh, 'Position', [0 0 .8 1]);
set(ppeh, 'backgroundColor', [1 1 1]);
set(ppeh, 'horizontalAlign', 'left');
set(ppeh, 'max', 2);
set(ppeh, 'min', 0);
set(ppeh, 'FontName', 'Courier', 'FontSize', editfontsize, 'FontUnits', editfontunits);
set(ppeh, 'string', code);

% add handle for the edit style to figure
setappdata(pph, 'superparent', parent); % superparent is the main ft_databrowser window
setappdata(pph, 'ppeh', ppeh);

end % function preproc_cfg1_cb


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function preproc_cfg2_cb(h, eventdata)
parent = get(h, 'parent');
superparent = getappdata(parent, 'superparent');
ppeh = getappdata(parent, 'ppeh');
code = get(ppeh, 'string');

% get rid of empty lines and white space
skip = [];
for iline = 1:numel(code)
  code{iline} = strtrim(code{iline});
  if isempty(code{iline})
    skip = [skip iline];
    continue
  end
  if code{iline}(1)=='%'
    skip = [skip iline];
    continue
  end
end
code(skip) = [];

if ~isempty(code)
  ispreproccfg = strncmp(code, 'cfg.preproc.',12);
  if ~all(ispreproccfg)
    errordlg('cfg-options must be specified as cfg.preproc.xxx', 'cfg.preproc editor', 'modal')
  end
  % eval the code
  for icomm = 1:numel(code)
    eval([code{icomm} ';']);
  end

  % check for cfg and output into the original appdata-window
  if ~exist('cfg', 'var')
    cfg = [];
    cfg.preproc = [];
  end
  maincfg = getappdata(superparent, 'cfg');
  maincfg.preproc = cfg.preproc;
  setappdata(superparent, 'cfg', maincfg)
end

close(parent)
redraw_cb(superparent)
end % function preproc_cfg2_cb


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enable_buttons(h, bool)
h = getparent(h);
buttons = findobj(h, 'style', 'pushbutton');
if bool
  set(buttons, 'enable', 'on');
else
  set(buttons, 'enable', 'off');
end
end % function enable_buttons


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enable_callbacks(h, bool)
h = getparent(h);
cfg = getappdata(h, 'cfg');
if bool
  set(h, 'ResizeFcn',             @resize_cb);
  set(h, 'KeyPressFcn',           @keypress_cb);
  set(h, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonDownFcn'});
  set(h, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonUpFcn'});
  set(h, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonMotionFcn'});
else
  set(h, 'ResizeFcn',             []);
  set(h, 'KeyPressFcn',           []);
  set(h, 'WindowButtonDownFcn',   []);
  set(h, 'WindowButtonUpFcn',     []);
  set(h, 'WindowButtonMotionFcn', []);
end
end % function enable_windowbuttoncallback


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function help_cb(h, eventdata)

fprintf('------------------------------------------------------------------------------------\n')
fprintf('You can use the following keyboard buttons in the databrowser\n')
fprintf('1-9                : select artifact type 1-9\n');
fprintf('shift 1-9          : select previous artifact of type 1-9 (does not work with numpad)\n');
fprintf('alt 1-9            : select next artifact of type 1-9\n');
fprintf('arrow-left         : previous trial\n');
fprintf('arrow-right        : next trial\n');
fprintf('shift arrow-up     : increase vertical scaling\n');
fprintf('shift arrow-down   : decrease vertical scaling\n');
fprintf('shift arrow-left   : increase horizontal scaling\n');
fprintf('shift arrow-down   : decrease horizontal scaling\n');
fprintf('t                  : open trial or segment selection dialog\n');
fprintf('c                  : open channel selection dialog\n');
fprintf('h                  : open horizontal scaling dialog\n');
fprintf('v                  : open vertical scaling dialog\n');
fprintf('p                  : open preproc editor\n');
fprintf('i                  : identify a specific channel\n');
fprintf('m                  : toggles between cfg.viewmode options\n');
fprintf('s                  : toggles between cfg.selectmode options\n');
fprintf('q                  : quit\n');
fprintf('------------------------------------------------------------------------------------\n')
fprintf('\n')
end % function help_cb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keypress_cb(h, eventdata)

% disable the buttons to prevent piling up of GUI events, they are enabled again at the end of the function
enable_buttons(h, false);

if (isempty(eventdata) && ft_platform_supports('matlabversion', -Inf, '2014a')) || isa(eventdata, 'matlab.ui.eventdata.ActionData')
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parsekeyboardevent(eventdata);
end

h = getparent(h);
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

switch key
  case {'1' '2' '3' '4' '5' '6' '7' '8' '9'}
    % switch to another artifact type
    opt.artsel = str2double(key);
    numart = size(opt.artdata.trial{1}, 1);
    if opt.artsel > numart
      fprintf('data has no artifact type %i \n', opt.artsel)
    else
      % bold the active one
      arth = findobj(h, 'tag', 'artifactui');
      arth = arth(end:-1:1); % order is reversed so reverse it again
      hsel = [1 2 3] + (opt.artsel-1) .*3 ;
      set(arth(hsel), 'fontweight', 'bold')
      % unbold the passive ones
      set(arth(setdiff(1:numel(arth),hsel)), 'fontweight', 'normal')
      % redraw
      setappdata(h, 'opt', opt);
      setappdata(h, 'cfg', cfg);
      fprintf('switching to the "%s" artifact\n', opt.artdata.label{opt.artsel});
      redraw_cb(h, eventdata);
    end
  case {'shift+1' 'shift+2' 'shift+3' 'shift+4' 'shift+5' 'shift+6' 'shift+7' 'shift+8' 'shift+9'}
    % go to previous artifact of the current type
    opt.artsel = str2double(key(end));
    numart = size(opt.artdata.trial{1}, 1);
    if opt.artsel > numart
      fprintf('data has no artifact type %i \n', opt.artsel)
    else
      % find the previous occuring artifact, keeping in mind that:
      % 1) artifacts can cross trial boundaries
      % 2) artifacts might not occur inside a trial boundary (when data is segmented differently than during artifact detection)
      % fetch trl representation of current artifact type
      arttrl = boolvec2trl(opt.artdata.trial{1}(opt.artsel,:));
      % discard artifacts in the future
      curvisend = opt.trlvis(opt.trlop,2);
      arttrl(arttrl(:,1) > curvisend,:) = [];
      % find nearest artifact by searching in each trl (we have to do this here everytime, because trlvis can change on the fly because of x-zooming)
      newtrlop = [];
      for itrlvis = opt.trlop-1:-1:1
        % is either the start or the end of any artifact present?
        if any(any(opt.trlvis(itrlvis,1)<=arttrl(:,1:2) & opt.trlvis(itrlvis,2)>=arttrl(:,1:2)))
          % if so, we're done
          newtrlop = itrlvis;
          break
        end
      end
      if isempty(newtrlop)
        fprintf('no earlier %s with "%s" artifact found\n', opt.trialviewtype, opt.artdata.label{opt.artsel});
      else
        fprintf('going to previous %s with "%s" artifact\n', opt.trialviewtype, opt.artdata.label{opt.artsel});
        opt.trlop = newtrlop;
        % other artifact type potentially selected, bold the active one
        arth = findobj(h, 'tag', 'artifactui');
        arth = arth(end:-1:1); % order is reversed so reverse it again
        hsel = [1 2 3] + (opt.artsel-1) .*3 ;
        set(arth(hsel), 'fontweight', 'bold')
        % unbold the passive ones
        set(arth(setdiff(1:numel(arth),hsel)), 'fontweight', 'normal')
        % export into fig and continue
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h, eventdata);
      end
    end
  case {'control+1' 'control+2' 'control+3' 'control+4' 'control+5' 'control+6' 'control+7' 'control+8' 'control+9' 'alt+1' 'alt+2' 'alt+3' 'alt+4' 'alt+5' 'alt+6' 'alt+7' 'alt+8' 'alt+9'}
    % go to next artifact
    opt.artsel = str2double(key(end));
    numart = size(opt.artdata.trial{1}, 1);
    if opt.artsel > numart
      fprintf('data has no artifact type %i \n', opt.artsel)
    else
      % find the next occuring artifact, keeping in mind that:
      % 1) artifacts can cross trial boundaries
      % 2) artifacts might not occur inside a trial boundary (when data is segmented differently than during artifact detection)
      % fetch trl representation of current artifact type
      arttrl = boolvec2trl(opt.artdata.trial{1}(opt.artsel,:));
      % discard artifacts in the past
      curvisbeg = opt.trlvis(opt.trlop,1);
      arttrl(arttrl(:,2) < curvisbeg,:) = [];
      % find nearest artifact by searching in each trl (we have to do this here everytime, because trlvis can change on the fly because of x-zooming)
      newtrlop = [];
      for itrlvis = opt.trlop+1:size(opt.trlvis,1)
        % is either the start or the end of any artifact present?
        if any(any(opt.trlvis(itrlvis,1)<=arttrl(:,1:2) & opt.trlvis(itrlvis,2)>=arttrl(:,1:2)))
          % if so, we're done
          newtrlop = itrlvis;
          break
        end
      end
      if isempty(newtrlop)
        fprintf('no later %s with "%s" artifact found\n', opt.trialviewtype, opt.artdata.label{opt.artsel});
      else
        fprintf('going to next %s with "%s" artifact\n', opt.trialviewtype, opt.artdata.label{opt.artsel});
        opt.trlop = newtrlop;
        % other artifact type potentially selected, bold the active one
        arth = findobj(h, 'tag', 'artifactui');
        arth = arth(end:-1:1); % order is reversed so reverse it again
        hsel = [1 2 3] + (opt.artsel-1) .*3 ;
        set(arth(hsel), 'fontweight', 'bold')
        % unbold the passive ones
        set(arth(setdiff(1:numel(arth),hsel)), 'fontweight', 'normal')
        % export into fig and continue
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h, eventdata);
      end
    end
  case 'leftarrow'
    opt.trlop = max(opt.trlop - 1, 1); % should not be smaller than 1
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    redraw_cb(h, eventdata);
  case 'rightarrow'
    opt.trlop = min(opt.trlop + 1, size(opt.trlvis,1)); % should not be larger than the number of trials
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    redraw_cb(h, eventdata);
  case 'uparrow'
    chansel = match_str(opt.hdr.label, cfg.channel);
    minchan = min(chansel);
    numchan = length(chansel);
    chansel = minchan - numchan : minchan - 1;
    if min(chansel)<1
      chansel = chansel - min(chansel) + 1;
    end
    % convert numeric array into cell-array with channel labels
    cfg.channel = opt.hdr.label(chansel);
    opt.changedchanflg = true; % trigger for redrawing channel labels and preparing layout again (see bug 2065 and 2878)
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    delete(findobj(h, 'tag', 'channellabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
    redraw_cb(h, eventdata);
  case 'downarrow'
    chansel = match_str(opt.hdr.label, cfg.channel);
    maxchan = max(chansel);
    numchan = length(chansel);
    chansel = maxchan + 1 : maxchan + numchan;
    if max(chansel)>length(opt.hdr.label)
      chansel = chansel - (max(chansel) - length(opt.hdr.label));
    end
    % convert numeric array into cell-array with channel labels
    cfg.channel = opt.hdr.label(chansel);
    opt.changedchanflg = true; % trigger for redrawing channel labels and preparing layout again (see bug 2065 and 2878)
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    delete(findobj(h, 'tag', 'channellabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
    redraw_cb(h, eventdata);
  case 'shift+leftarrow'
    cfg.blocksize = cfg.blocksize*sqrt(2);
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    definetrial_cb(h, eventdata);
    redraw_cb(h, eventdata);
  case 'shift+rightarrow'
    cfg.blocksize = cfg.blocksize/sqrt(2);
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    definetrial_cb(h, eventdata);
    redraw_cb(h, eventdata);
  case 'shift+uparrow'
    cfg.ylim = cfg.ylim/sqrt(2);
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    redraw_cb(h, eventdata);
  case 'shift+downarrow'
    cfg.ylim = cfg.ylim*sqrt(2);
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    redraw_cb(h, eventdata);
  case 'p'
    preproc_cfg1_cb(h, eventdata);
  case 'q'
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    cleanup_cb(h);
  case 't'
    % select the trial to display
    if ~strcmp(opt.trialviewtype, 'trialsegment')
      str = sprintf('%s to display (current trial = %d/%d)', opt.trialviewtype, opt.trlop, size(opt.trlvis,1));
    else
      str = sprintf('segment to display (current segment = %d/%d)', opt.trlop, size(opt.trlvis,1));
    end
    response = inputdlg(str, 'specify', 1, {num2str(opt.trlop)});
    if ~isempty(response)
      opt.trlop = str2double(response);
      opt.trlop = min(opt.trlop, size(opt.trlvis,1)); % should not be larger than the number of trials
      opt.trlop = max(opt.trlop, 1); % should not be smaller than 1
      setappdata(h, 'opt', opt);
      setappdata(h, 'cfg', cfg);
      redraw_cb(h, eventdata);
    end
  case 'h'
    % select the horizontal scaling
    response = inputdlg('horizontal scale', 'specify', 1, {num2str(cfg.blocksize)});
    if ~isempty(response)
      cfg.blocksize = str2double(response);
      setappdata(h, 'opt', opt);
      setappdata(h, 'cfg', cfg);
      definetrial_cb(h, eventdata);
      redraw_cb(h, eventdata);
    end
  case 'v'
    % select the vertical scaling
    response = inputdlg('vertical scale, [ymin ymax], ''maxabs'' or ''maxmin''', 'specify', 1, {['[ ' num2str(cfg.ylim) ' ]']});
    if ~isempty(response)
      if strcmp(response, 'maxmin')
        minval = min(opt.curdata.trial{1}(:));
        maxval = max(opt.curdata.trial{1}(:));
        cfg.ylim = [minval maxval];
      elseif strcmp(response, 'maxabs')
        minval = min(opt.curdata.trial{1}(:));
        maxval = max(opt.curdata.trial{1}(:));
        cfg.ylim = [-max(abs([minval maxval])) max(abs([minval maxval]))];
      else
        % convert to string and add brackets, just to ensure that str2num will work
        tmp = str2num(['[' response{1} ']']);
        if numel(tmp)==2
          cfg.ylim = tmp;
        else
          ft_warning('incorrect specification of cfg.ylim, not changing the limits for the vertical axes')
        end
      end
      setappdata(h, 'opt', opt);
      setappdata(h, 'cfg', cfg);
      redraw_cb(h, eventdata);
    end
  case 'c'
    % select channels
    select = match_str(opt.hdr.label, cfg.channel);
    select = select_channel_list(opt.hdr.label, select);
    cfg.channel = opt.hdr.label(select);
    opt.changedchanflg = true; % trigger for redrawing channel labels and preparing layout again (see bug 2065 and 2878)
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    delete(findobj(h, 'tag', 'channellabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
    redraw_cb(h, eventdata);
  case 'i'
    % click in the figure to highlight and get the name of the nearest channel
    delete(findobj(h, 'tag', 'identifiedchannel'));
    fprintf('click in the figure to identify the name of the closest channel\n');
    val = ginput(1);
    xpos = val(1);
    if strcmp(cfg.viewmode, 'butterfly')
      if val>0.5
        ypos = 0.1;
      else
        ypos = 0.9;
      end
      % transform 'val' to match data
      val(1) = val(1) * range(opt.hlim) + opt.hlim(1);
      val(2) = val(2) * range(opt.vlim) + opt.vlim(1);
      label  = val2nearestchan(opt.curdata, val);
      datindx    = match_str(opt.curdata.label, label);
      layoutindx = 1; % the butterfly layout only contains a single dummy channel that is used for all
    elseif strcmp(cfg.viewmode, 'vertical') || strcmp(cfg.viewmode, 'component')
      % find channel identity by extracting timecourse objects and finding the time course closest to the cursor
      % this is a lot easier than the reverse, determining the y value of each time course scaled by the layout and vlim
      tcobj   = findobj(h, 'tag', 'timecourse');
      tmpydat = get(tcobj, 'ydata');
      tmpydat = cat(1,tmpydat{:});
      tmpydat = tmpydat(end:-1:1,:); % order of timecourse objects is reverse of channel order
      tmpxdat = get(tcobj(1), 'xdata');
      % first find the closest sample number along the horizontal direction
      xsmp = nearest(tmpxdat, val(1));
      % then find the closest value over all channels, this gives the channel number
      datindx    = nearest(tmpydat(:,xsmp), val(2));
      label      = opt.curdata.label{datindx};
      % decide whether to place the label below or above the line
      layoutindx = match_str(opt.layouttime.label, label);
      if val(2)>opt.layouttime.pos(layoutindx,2)
        ypos = opt.layouttime.pos(layoutindx,2) - opt.layouttime.height(layoutindx)/2;
      else
        ypos = opt.layouttime.pos(layoutindx,2) + opt.layouttime.height(layoutindx)/2;
      end
    end
    fprintf('identified channel name: %s\n',label);
    redraw_cb(h, eventdata);

    ft_plot_text(xpos, ypos, label, 'FontSize', cfg.fontsize, 'FontUnits', cfg.fontunits, 'tag', 'identifiedchannel', 'FontSize', cfg.fontsize, 'FontUnits', cfg.fontunits);
    if ~ishold
      hold on
      ft_plot_vector(opt.curdata.time{1}, opt.curdata.trial{1}(datindx,:)', 'box', false, 'tag', 'identifiedchannel', 'hpos', opt.layouttime.pos(layoutindx,1), 'vpos', opt.layouttime.pos(layoutindx,2), 'width', opt.layouttime.width(layoutindx), 'height', opt.layouttime.height(layoutindx), 'hlim', opt.hlim, 'vlim', opt.vlim, 'color', 'k', 'linewidth', 2);
      hold off
    else
      ft_plot_vector(opt.curdata.time{1}, opt.curdata.trial{1}(datindx,:)', 'box', false, 'tag', 'identifiedchannel', 'hpos', opt.layouttime.pos(layoutindx,1), 'vpos', opt.layouttime.pos(layoutindx,2), 'width', opt.layouttime.width(layoutindx), 'height', opt.layouttime.height(layoutindx), 'hlim', opt.hlim, 'vlim', opt.vlim, 'color', 'k', 'linewidth', 2);
    end
  case 's'
    % toggle between selectmode options
    curstate = find(strcmp(cfg.selectmode, {'markartifact', 'markpeakevent', 'marktroughevent'}));
    if curstate == 1
      cfg.selectmode = 'markpeakevent'; % go to the next
    elseif curstate == 2
      cfg.selectmode = 'marktroughevent'; % go to the next
    elseif curstate == 3
      cfg.selectmode = 'markartifact'; % go to the first
    end
    fprintf('switching to selectmode = %s\n', cfg.selectmode);
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    redraw_cb(h, eventdata);
  case 'm'
    % it is OK to have some overlap between the traces in vertical and component viewmodes
    % but butterfly plots should stay within the min/max boundaries
    scale_adjust = 5;
    switch cfg.viewmode
      case 'butterfly'
        % toggle between viewmodes
        cfg.viewmode = 'vertical';
        cfg.ylim = cfg.ylim / scale_adjust;
        opt.changedchanflg = true;
        delete(findobj(h, 'tag', 'channellabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h, eventdata);
      case 'vertical'
        % toggle between viewmodes
        cfg.viewmode = 'butterfly';
        cfg.ylim = cfg.ylim * scale_adjust;
        opt.changedchanflg = true;
        delete(findobj(h, 'tag', 'channellabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h, eventdata);
      otherwise
        % do nothing
    end % switch
  case 'control+control'
    % do nothing
  case 'shift+shift'
    % do nothing
  case 'alt+alt'
    % do nothing
  otherwise
    help_cb(h, eventdata);
end

enable_buttons(h, true);
uiresume(h);
end % function keypress_cb


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resize_cb(h, eventdata)
h = getparent(h);
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

if strcmp(cfg.viewmode, 'butterfly')
  % no reason to redraw
else
  opt.changedchanflg = true;
  setappdata(h, 'opt', opt);
  redraw_cb(h, eventdata);
end

end % function resize_cb


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function redraw_cb(h, eventdata)

% disable the interactive callbacks to prevent piling up of GUI events, they are enabled again at the end of the function
enable_callbacks(h, false);

h = getparent(h);
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

% ensure that the calling figure is in the front
if h~=gcf
  figure(h);
end

delete(findobj(h, 'tag', 'timecourse'));
delete(findobj(h, 'tag', 'identifiedchannel'));
delete(findobj(h, 'tag', 'selectedrange'));
% not removing channel labels, they cause the bulk of redrawing time for the slow text function (note, interpreter = none hardly helps), see bug 2065
% delete(findobj(h, 'tag', 'channellabel'));



% temporarily store the originally selected list of channels
userchan    = cfg.channel;

% add the clamped channels at the end of the list
cfg.channel = setdiff(cfg.channel, cfg.channelclamped, 'stable');
cfg.channel = vertcat(cfg.channel(:), cfg.channelclamped(:));

% if the number of channels changes, includes past channels
if numel(cfg.channel)<numel(userchan) + numel(cfg.channelclamped)
  addchan = numel(userchan) + numel(cfg.channelclamped) - numel(cfg.channel); % Gets the number of channels to add.
  lastchan = min(match_str(opt.hdr.label, cfg.channel)) - 1; % Gets the previous channel.

  % Adds up to n more channels.
  cfg.channel = union(opt.hdr.label(max(lastchan - addchan + 1, 1): lastchan), cfg.channel);
end

begsample = opt.trlvis(opt.trlop, 1);
endsample = opt.trlvis(opt.trlop, 2);
offset    = opt.trlvis(opt.trlop, 3);
chanindx  = match_str(opt.hdr.label, cfg.channel);

if isempty(opt.orgdata)
  dat = ft_read_data(cfg.datafile, 'header', opt.hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', ~istrue(cfg.continuous), 'dataformat', cfg.dataformat, opt.headeropt{:});
else
  dat = ft_fetch_data(opt.orgdata, 'header', opt.hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'allowoverlap', cfg.allowoverlap, 'skipcheckdata', true);
end

% fetch only the artifacts in the current time window, they are represented as booleans, with one channel per artifact type
artdat = ft_fetch_data(opt.artdata, 'begsample', begsample, 'endsample', endsample);
artlab = opt.artdata.label;

if ~isempty(opt.event) && isstruct(opt.event)
  % select those events that overlap with the current time window
  event    = opt.event;
  begevent = [event(:).sample];
  endevent = [event(:).sample] + [event(:).duration];
  event    = event(begevent<=endsample & endevent>=begsample);
else
  event = [];
end

% convert the data to another numeric precision, i.e. double, single or int32
if ~isempty(cfg.precision)
  dat = cast(dat, cfg.precision);
end

% apply preprocessing and determine the time axis
tim = offset2time(offset, opt.fsample, size(dat,2));
finitevals = isfinite(sum(dat));
if all(~finitevals)
  lab = opt.hdr.label(chanindx);
elseif any(~finitevals)
  % loop across the chunks to avoid nan-related warnings from preproc and potential loss of data
  begsmp = find( ~[0 finitevals] & [finitevals 0] );
  endsmp = find( ~[finitevals 0] & [0 finitevals] ) - 1;
  for k = 1:numel(begsmp)
    [dat(:,begsmp(k):endsmp(k)), lab] = preproc(dat(:,begsmp(k):endsmp(k)), opt.hdr.label(chanindx), tim(begsmp(k):endsmp(k)), cfg.preproc);
  end
else
  [dat, lab, tim] = preproc(dat, opt.hdr.label(chanindx), tim, cfg.preproc);
end

% add NaNs to data for plotting purposes. NaNs will be added when requested horizontal scaling is longer than the data.
nsamplepad = opt.nanpaddata(opt.trlop);
if nsamplepad>0
  dat = [dat NaN(numel(lab), opt.nanpaddata(opt.trlop))];
  tim = [tim linspace(tim(end),tim(end)+nsamplepad*mean(diff(tim)),nsamplepad)];  % possible machine precision error here
end

% make a single-trial data structure for the current data
opt.curdata.label      = lab;
opt.curdata.time{1}    = tim;
opt.curdata.trial{1}   = dat;
opt.curdata.hdr        = opt.hdr;
opt.curdata.fsample    = opt.fsample;
opt.curdata.sampleinfo = [begsample endsample];
opt.curdata = copyfields(opt.orgdata, opt.curdata, {'grad', 'elec', 'opto'});
% remove the local copy of the data fields
clear lab tim dat

fn = fieldnames(cfg);
tmpcfg = keepfields(cfg, fn(endsWith(fn, 'scale') | startsWith(fn, 'mychan') | strcmp(fn, 'channel')));
tmpcfg.parameter = 'trial';
opt.curdata = chanscale_common(tmpcfg, opt.curdata);

% make a local copy (again) of the data fields
lab = opt.curdata.label;
tim = opt.curdata.time{1};
dat = opt.curdata.trial{1};

% this is a trigger for redrawing channel labels and preparing layout again (see bug 2065 and 2878)
changedchanflg = opt.changedchanflg;
opt.changedchanflg = false;

if strcmp(cfg.viewmode, 'butterfly')
  % the timecourse layout is always the same
  layouttime = [];
  layouttime.label = {'dummy'};
  layouttime.pos = [0.5 0.5];
  layouttime.width = 1;
  layouttime.height = 1;
  opt.layouttime = layouttime;
else
  % the timecourse layout needs to be reconstructed whenever the channel selection changes
  if changedchanflg % trigger for redrawing channel labels and preparing layout again (see bug 2065 and 2878)
    tmpcfg = keepfields(cfg, {'channel', 'columns', 'rows', 'commentpos', 'scalepos', 'elec', 'grad', 'opto', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
    tmpcfg.layout = 'vertical';
    tmpcfg.skipcomnt = 'yes';
    tmpcfg.skipscale = 'yes';
    tmpcfg.showcallinfo = 'no';
    tmpcfg.trackcallinfo = 'no';
    if ~isempty(opt.orgdata)
      tmpdata = ft_checkdata(opt.orgdata, 'datatype', 'raw'); % remove the component information from the data
      opt.layouttime = ft_prepare_layout(tmpcfg, tmpdata);
    else
      opt.layouttime = ft_prepare_layout(tmpcfg);
    end
  end
end

% determine the position of the channel/component labels relative to the real axes
if strcmp(cfg.viewmode, 'component')
  labelx = opt.layouttime.pos(:,1) - opt.layouttime.width/2 - 2*opt.layouttime.height;
else
  labelx = opt.layouttime.pos(:,1) - opt.layouttime.width/2 - 0.01;
end
labely = opt.layouttime.pos(:,2);

% determine the total extent of all virtual axes relative to the real axes
ax(1) = min(opt.layouttime.pos(:,1) - opt.layouttime.width/2);
ax(2) = max(opt.layouttime.pos(:,1) + opt.layouttime.width/2);
ax(3) = min(opt.layouttime.pos(:,2) - opt.layouttime.height/2);
ax(4) = max(opt.layouttime.pos(:,2) + opt.layouttime.height/2);
% add white space to bottom and top so channels are not out-of-axis for the majority
% NOTE there is a second spot where this is done below, specifically for viewmode = component (also need to be here), which should be kept the same as this
if strcmp(cfg.viewmode, 'vertical') || strcmp(cfg.viewmode, 'component')
  % determine amount of vertical padding using cfg.verticalpadding
  if ~isnumeric(cfg.verticalpadding) && strcmp(cfg.verticalpadding, 'auto')
    % determine amount of padding using the number of channels
    if numel(cfg.channel)<=6
      wsfac = 0;
    elseif numel(cfg.channel)>6 && numel(cfg.channel)<=10
      wsfac = 0.01 *  (ax(4)-ax(3));
    else
      wsfac = 0.02 *  (ax(4)-ax(3));
    end
  else
    wsfac = cfg.verticalpadding * (ax(4)-ax(3));
  end
  ax(3) = ax(3) - wsfac;
  ax(4) = ax(4) + wsfac;
end
axis(ax)

% determine a single local axis that encompasses all channels
% this is in relative figure units
opt.hpos   = (ax(1)+ax(2))/2;
opt.vpos   = (ax(3)+ax(4))/2;
opt.width  = ax(2)-ax(1);
opt.height = ax(4)-ax(3);

% these determine the scaling inside the virtual axes
% the hlim will be in seconds, the vlim will be in Tesla or Volt
opt.hlim = [tim(1) tim(end)];
opt.vlim = cfg.ylim;

if strcmp(cfg.plotartifacts, 'yes')
  delete(findobj(h, 'tag', 'artifact'));

  % ensure that the currently active artifact is plotted the last (and hence on top)
  ordervec = 1:length(artlab);
  ordervec(opt.artsel) = [];
  ordervec(end+1) = opt.artsel;

  % save stuff to be able to shift the labels downwards when they occur close to each other
  artifacttime = NaN(1,1);
  artifactshift = zeros(1,1);

  i = 0;
  for j = ordervec
    tmp = diff([0 artdat(j,:) 0]);
    artbeg = find(tmp==+1);
    artend = find(tmp==-1) - 1;

    % construct the text label that will be shown with the artifacts
    switch cfg.plotartifactlabels
      case 'type'
        artifactlabel = artlab{j};
      otherwise
        artifactlabel = '';
    end

    for k=1:numel(artbeg)
      i = i + 1; % it is not a simple loop, there are multiple types of artifacts, times multiple occurences
      artifacttime(i) = tim(artbeg(k)) + opt.hlim(1);
      xpos = [tim(artbeg(k)) tim(artend(k))];

      if xpos(1)==xpos(2)
        % plot it as a line when it has no duration
        lh = ft_plot_line(xpos, [-1 1], 'tag', 'artifact', 'hpos', opt.hpos, 'vpos', opt.vpos, 'width', opt.width, 'height', opt.height, 'hlim', opt.hlim, 'vlim', [-1 1], 'color', opt.artifactcolors(j,:));
      else
        % plot it as a box when it has a duration, pad it on either side with half a sample
        xpos = xpos + [-.5 +.5]./opt.fsample;
        lh = ft_plot_box([xpos -1 1], 'tag', 'artifact', 'hpos', opt.hpos, 'vpos', opt.vpos, 'width', opt.width, 'height', opt.height, 'hlim', opt.hlim, 'vlim', [-1 1], 'edgecolor', 'none', 'facecolor', opt.artifactcolors(j,:), 'facealpha', cfg.artifactalpha);
      end

      % store this data in the line object so that it can be displayed with cb_datacursortext
      setappdata(lh, 'datacursor_linetype', 'artifact');
      setappdata(lh, 'datacursor_artifacttime', artifacttime(i));
      setappdata(lh, 'datacursor_artifactlabel', artlab{j}); % always use the artifact type in the datacursor

      % check for events that are within 1/10th of the horizontal time axis and give the event label some vertical shift
      closeevent = find(artifacttime(i)>(artifacttime(1:i-1)-diff(opt.hlim)/10) & artifacttime(i)<(artifacttime(1:i-1)+diff(opt.hlim)/10));
      emptyshift = find(diff(sort(artifactshift(closeevent)))>1,1); % find a shift that has not been used by other close events
      if ~isempty(emptyshift)
        artifactshift(i) = emptyshift;
      elseif ~any(artifactshift(closeevent)==0) % restart at 0 when no other close event is plotted at the highest level
        artifactshift(i) = 0;
      elseif closeevent
        artifactshift(i) = max(artifactshift(closeevent))+1;
      else
        artifactshift(i) = 0;
      end

      % plot the artifact label
      ft_plot_text(artifacttime(i), -1+artifactshift(i)*.05, artifactlabel, 'tag', 'artifact', 'color',  opt.artifactcolors(j,:), 'hpos', opt.hpos, 'vpos', opt.vpos, 'width', opt.width, 'height', opt.height, 'hlim', opt.hlim, 'vlim', [-1 1], 'FontSize', cfg.fontsize, 'FontUnits', cfg.fontunits, 'horizontalalignment', 'left', 'verticalalignment', 'bottom');

    end % for each artifact

  end % for each of the artifact channels
end % if plot artifacts


if strcmp(cfg.plotevents, 'yes')
  delete(findobj(h, 'tag', 'event'));

  % save stuff to be able to shift the labels downwards when they occur close to each other
  eventtime  = NaN(1,numel(event));
  eventshift = zeros(1,numel(event));

  % plot a line or a box for each of the events
  for i = 1:numel(event)
    eventtype     = event(i).type;
    eventvalue    = event(i).value;
    eventduration = event(i).duration;
    eventoffset   = event(i).offset;

    % construct the text label that will be shown with the event
    switch cfg.ploteventlabels
      case 'type=value'
        if ~isempty(eventvalue)
          eventlabel = sprintf('%s=%s', eventtype, num2str(eventvalue)); % value can be both number and string
        else
          eventlabel = eventtype;
        end
      case 'type'
        eventlabel = eventtype;
      case 'value'
        if ~isempty(eventvalue)
          eventlabel = sprintf('%s', num2str(eventvalue)); % value can be both number and string
        else
          eventlabel = '';
        end
      otherwise
        ft_warning('unsupported specification of cfg.ploteventlabels');
        eventlabel = '';
    end

    if isempty(opt.eventcolors)
      eventcol = 'k';
    else
      eventcol = opt.eventcolors(strcmp(opt.eventtypes, eventtype),:);
    end

    % compute the time that the event starts
    eventtime(i) = (event(i).sample-begsample)/opt.fsample + opt.hlim(1);
    xpos = [eventtime(i) eventtime(i)+eventduration/opt.fsample];

    if length(xpos)>1 && xpos(1)==xpos(2)
      % plot it as a line when it has no duration
      lh = ft_plot_line(xpos, [-1 1], 'tag', 'event', 'hpos', opt.hpos, 'vpos', opt.vpos, 'width', opt.width, 'height', opt.height, 'hlim', opt.hlim, 'vlim', [-1 1], 'color', eventcol);
    else
      % plot it as a box when it has a duration, pad it on either side with half a sample
      xpos = xpos + [-.5 +.5]./opt.fsample;
      lh = ft_plot_box([xpos -1 1], 'tag', 'event', 'hpos', opt.hpos, 'vpos', opt.vpos, 'width', opt.width, 'height', opt.height, 'hlim', opt.hlim, 'vlim', [-1 1], 'edgecolor', eventcol, 'facecolor', eventcol, 'facealpha', cfg.eventalpha);
      if eventoffset<0 && eventoffset>-eventduration
        % add a vertical line, but only if it falls inside the box
        xpos = [eventtime(i) eventtime(i)] - eventoffset/opt.fsample;
        ft_plot_line(xpos, [-1 1], 'tag', 'event', 'hpos', opt.hpos, 'vpos', opt.vpos, 'width', opt.width, 'height', opt.height, 'hlim', opt.hlim, 'vlim', [-1 1], 'color', eventcol);
      end
    end

    % store this data in the line object so that it can be displayed with cb_datacursortext
    setappdata(lh, 'datacursor_linetype', 'event');
    setappdata(lh, 'datacursor_eventtime', eventtime(i));
    setappdata(lh, 'datacursor_eventlabel', sprintf('%s=%s', eventtype, num2str(eventvalue))); % always use type=value in the datacursor

    % check for events that are within 1/10th of the horizontal time axis and give the event label some vertical shift
    closeevent = find(eventtime(i)>(eventtime(1:i-1)-diff(opt.hlim)/10) & eventtime(i)<(eventtime(1:i-1)+diff(opt.hlim)/10));
    emptyshift = find(diff(sort(eventshift(closeevent)))>1,1); % find a shift that has not been used by other close events
    if ~isempty(emptyshift)
      eventshift(i) = emptyshift;
    elseif ~any(eventshift(closeevent)==0) % restart at 0 when no other close event is plotted at the highest level
      eventshift(i) = 0;
    elseif closeevent
      eventshift(i) = max(eventshift(closeevent))+1;
    else
      eventshift(i) = 0;
    end

    % plot the event label
    ft_plot_text(eventtime(i), 1-eventshift(i)*.05, eventlabel, 'tag', 'event', 'color', eventcol, 'hpos', opt.hpos, 'vpos', opt.vpos, 'width', opt.width, 'height', opt.height, 'hlim', opt.hlim, 'vlim', [-1 1], 'FontSize', cfg.fontsize, 'FontUnits', cfg.fontunits, 'horizontalalignment', 'left', 'verticalalignment', 'top');
  end % for numel(event)

end % if plot events

if strcmp(cfg.viewmode, 'butterfly')
  ft_plot_vector(tim, dat', 'box', false, 'tag', 'timecourse', 'hpos', opt.layouttime.pos(1,1), 'vpos', opt.layouttime.pos(1,2), 'width', opt.layouttime.width(1), 'height', opt.layouttime.height(1), 'hlim', opt.hlim, 'vlim', opt.vlim, 'color', opt.linecolor(chanindx,:), 'linewidth', cfg.linewidth, 'style', cfg.linestyle);
  set(gca, 'FontSize', cfg.axisfontsize, 'FontUnits', cfg.axisfontunits);

  % two ticks per channel
  yTick = sort([
    opt.layouttime.pos(:,2)+(opt.layouttime.height/2)
    opt.layouttime.pos(:,2)+(opt.layouttime.height/4)
    opt.layouttime.pos(:,2)
    opt.layouttime.pos(:,2)-(opt.layouttime.height/4)
    opt.layouttime.pos(:,2)-(opt.layouttime.height/2)
    ]); % sort

  yTickLabel = {num2str(yTick.*range(opt.vlim) + opt.vlim(1))};

  set(gca, 'yTick', yTick, 'yTickLabel', yTickLabel)

elseif strcmp(cfg.viewmode, 'vertical') || strcmp(cfg.viewmode, 'component')
  % determine channel indices into data outside of loop
  laysels = match_str(opt.layouttime.label, opt.hdr.label);

  % delete old chan labels before renewing, if they need to be renewed
  if changedchanflg % trigger for redrawing channel labels and preparing layout again (see bug 2065 and 2878)
    delete(findobj(h, 'tag', 'channellabel'));
  end

  % determine how many labels to skip in case of 'some'
  if strcmp(cfg.plotlabels, 'some') && strcmp(cfg.fontunits, 'points')
    % determine number of labels to plot by estimating overlap using current figure size
    % the idea is that figure height in pixels roughly corresponds to the amount of letters at cfg.fontsize (points) put above each other without overlap
    figheight = get(h, 'Position');
    figheight = figheight(4);
    labheight = cfg.fontsize+6; % 6 added, so that labels are not too close together (i.e. overlap if font was 6 points bigger)
    labskip = ceil(numel(chanindx)*labheight/figheight);
  elseif strcmp(cfg.plotlabels, 'some')
    labskip = 10;
  else
    % do not skip any labels
    labskip = 1;
  end

  for i = 1:length(chanindx)
    datsel = i; % this is not the original channel number, but the one from the selection
    laysel = laysels(i);

    if ~isempty(datsel) && ~isempty(laysel)
      % only plot chanlabels when necessary
      if changedchanflg % trigger for redrawing channel labels and preparing layout again (see bug 2065 and 2878)

        if strcmp(cfg.plotlabels, 'yes') || (strcmp(cfg.plotlabels, 'some') && mod(i,labskip)==0)
          ft_plot_text(labelx(laysel), labely(laysel), opt.hdr.label(chanindx(i)), 'tag', 'channellabel', 'HorizontalAlignment', 'right', 'FontSize', cfg.fontsize, 'FontUnits', cfg.fontunits);
          set(gca, 'FontSize', cfg.axisfontsize, 'FontUnits', cfg.axisfontunits);
        end
      end

      lh = ft_plot_vector(tim, dat(datsel, :)', 'box', false, 'tag', 'timecourse', 'hpos', opt.layouttime.pos(laysel,1), 'vpos', opt.layouttime.pos(laysel,2), 'width', opt.layouttime.width(laysel), 'height', opt.layouttime.height(laysel), 'hlim', opt.hlim, 'vlim', opt.vlim, 'color', opt.linecolor(chanindx(i),:), 'linewidth', cfg.linewidth, 'style', cfg.linestyle);

      % store this data in the line object so that it can be displayed with cb_datacursortext
      setappdata(lh, 'datacursor_linetype', 'channel');
      setappdata(lh, 'datacursor_label', opt.hdr.label(chanindx(i)));
      setappdata(lh, 'datacursor_xdata', tim);
      setappdata(lh, 'datacursor_ydata', dat(datsel,:));
    end
  end

  % plot yticks
  if length(chanindx)>6
    % plot yticks at each label in case adaptive labeling is used (cfg.plotlabels = 'some')
    % otherwise, use the old ytick plotting based on hard-coded number of channels
    if strcmp(cfg.plotlabels, 'some')
      yTick = sort(labely(mod(chanindx,labskip)==0), 'ascend'); % sort is required, yticks should be increasing in value
      yTickLabel = [];
    else
      if length(chanindx)>19
        % no space for yticks
        yTick = [];
        yTickLabel = [];
      elseif length(chanindx)> 6
        % one tick per channel
        yTick = sort([
          opt.layouttime.pos(:,2)+(opt.layouttime.height(laysel)/4)
          opt.layouttime.pos(:,2)-(opt.layouttime.height(laysel)/4)
          ]);
        yTickLabel = {[.25 .75] .* range(opt.vlim) + opt.vlim(1)};
      end
    end
  else
    % two ticks per channel
    yTick = sort([
      opt.layouttime.pos(:,2)+(opt.layouttime.height(laysel)/2)
      opt.layouttime.pos(:,2)+(opt.layouttime.height(laysel)/4)
      opt.layouttime.pos(:,2)-(opt.layouttime.height(laysel)/4)
      opt.layouttime.pos(:,2)-(opt.layouttime.height(laysel)/2)
      ]); % sort
    yTickLabel = {[.0 .25 .75 1] .* range(opt.vlim) + opt.vlim(1)};
  end
  yTickLabel = repmat(yTickLabel, 1, length(chanindx));
  set(gca, 'yTick', yTick, 'yTickLabel', yTickLabel);

else
  % the following is implemented for other viewmodes
  % for example when the channels are organized according to a CTF151 layout

  % determine channel indices into data outside of loop
  laysels = match_str(opt.layouttime.label, opt.hdr.label);

  for i = 1:length(chanindx)
    datsel = i;
    laysel = laysels(i);

    if ~isempty(datsel) && ~isempty(laysel)
      lh = ft_plot_vector(tim, dat(datsel, :)', 'box', false, 'tag', 'timecourse', 'hpos', opt.layouttime.pos(laysel,1), 'vpos', opt.layouttime.pos(laysel,2), 'width', opt.layouttime.width(laysel), 'height', opt.layouttime.height(laysel), 'hlim', opt.hlim, 'vlim', opt.vlim, 'color', opt.linecolor(chanindx(i),:), 'linewidth', cfg.linewidth, 'style', cfg.linestyle);

      % store this data in the line object so that it can be displayed with cb_datacursortext
      setappdata(lh, 'datacursor_linetype', 'channel');
      setappdata(lh, 'datacursor_label', opt.hdr.label(chanindx(i)));
      setappdata(lh, 'datacursor_xdata', tim);
      setappdata(lh, 'datacursor_ydata', dat(datsel,:));
    end
  end

  % ticks are not supported with such a layout
  yTick = [];
  yTickLabel = [];
  yTickLabel = repmat(yTickLabel, 1, length(chanindx));
  set(gca, 'yTick', yTick, 'yTickLabel', yTickLabel);

end % if strcmp viewmode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add the ticks along the horizontal axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(strcmp(cfg.viewmode, {'butterfly', 'component', 'vertical'}))
  nticks = 11;
  xTickLabel = cellstr(num2str( linspace(tim(1), tim(end), nticks)' , '%1.2f'))';
  xTick = linspace(ax(1), ax(2), nticks);
  if nsamplepad>0
    nlabindat = sum(linspace(tim(1), tim(end), nticks) < tim(end-nsamplepad));
    xTickLabel(nlabindat+1:end) = repmat({' '}, [1 nticks-nlabindat]);
  end
  set(gca, 'xTick', xTick, 'xTickLabel', xTickLabel)
  xlabel('time');
else
  set(gca, 'xTick', [], 'xTickLabel', [])
end

if strcmp(cfg.viewmode, 'component')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % add the component topographies
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % determine the layout for each of the component topographies
  layoutcomp.pos(:,1)  = opt.layouttime.pos(:,1) - opt.layouttime.width/2 - opt.layouttime.height;
  layoutcomp.pos(:,2)  = opt.layouttime.pos(:,2);
  layoutcomp.width     = opt.layouttime.height;
  layoutcomp.height    = opt.layouttime.height;
  layoutcomp.label     = opt.layouttime.label;

  if ~isequal(opt.chanindx, chanindx)
    opt.chanindx = chanindx;

    delete(findobj(h, 'tag', 'topography'));

    % determine the position of each of the original channels for the topgraphy
    [sel1, sel2] = match_str(opt.orgdata.topolabel, opt.layouttopo.label);
    chanx = opt.layouttopo.pos(sel2,1);
    chany = opt.layouttopo.pos(sel2,2);

    if strcmp(cfg.compscale, 'global')
      % compute scaling factors for all components that are plotted together
      zmin = nan(1, length(chanindx));
      zmax = nan(1, length(chanindx));
      for i=1:length(chanindx)
        zmin(i) = min(opt.orgdata.topo(sel1,chanindx(i)));
        zmax(i) = max(opt.orgdata.topo(sel1,chanindx(i)));
      end

      if strcmp(cfg.zlim, 'maxmin')
        zmin = min(zmin);
        zmax = max(zmax);
      elseif strcmp(cfg.zlim, 'maxabs')
        zmax = max(abs([zmin zmax]));
        zmin = -zmax;
      elseif strcmp(cfg.zlim, 'zeromax')
        zmin = 0;
        zmax = max(zmax);
      elseif strcmp(cfg.zlim, 'minzero')
        zmin = min(zmin);
        zmax = 0;
      elseif isnumeric(cfg.zlim) && numel(cfg.zlim)==2
        zmin = cfg.zlim(1);
        zmax = cfg.zlim(2);
      else
        ft_error('unsupported specification of cfg.zlim');
      end
    end % if cfg.compscale

    for i=1:length(chanindx)
      % plot the topography of this component
      laysel = match_str(opt.layouttime.label, opt.hdr.label(chanindx(i)));
      chanz = opt.orgdata.topo(sel1,chanindx(i));
      if all(chanz==0)
        % this is most likely a channel that is not a component, and has been added later on
        continue;
      end

      if strcmp(cfg.compscale, 'local')
        % compute scaling factors for each individual component
        if strcmp(cfg.zlim, 'maxmin')
          zmin = min(chanz);
          zmax = max(chanz);
        elseif strcmp(cfg.zlim, 'maxabs')
          zmax = max(abs(chanz));
          zmin = -zmax;
        elseif strcmp(cfg.zlim, 'zeromax')
          zmin = 0;
          zmax = max(chanz);
        elseif strcmp(cfg.zlim, 'minzero')
          zmin = min(chanz);
          zmax = 0;
        elseif isnumeric(cfg.zlim) && numel(cfg.zlim)==2
          zmin = cfg.zlim(1);
          zmax = cfg.zlim(2);
        else
          ft_error('unsupported specification of cfg.zlim');
        end
      end % if cfg.compscale

      % scaling
      chanz = (chanz - zmin) ./  (zmax- zmin);

      % layouttopo is the actual 2D topographic channel layout, in pixel units for .mat files
      % layoutcomp is a vertical layout determining where to plot each component topography, with one entry per component

      plotopt = {
        'interpmethod', cfg.interpolation, ...
        'interplim',    cfg.interplimits, ...
        'gridscale',    cfg.gridscale, ...
        'outline',      opt.layouttopo.outline, ...
        'mask',         opt.layouttopo.mask, ...
        'shading',      cfg.shading, ...
        'isolines',     cfg.contournum, ...
        'tag',          'topography', ...
        'hpos',         layoutcomp.pos(laysel,1), ...
        'vpos',         layoutcomp.pos(laysel,2), ...
        'width',        layoutcomp.width(laysel), ...
        'height',       layoutcomp.height(laysel)};

      ft_plot_topo(chanx, chany, chanz, plotopt{:});
      % axis equal
      % drawnow
    end % for chanindx

    caxis([0 1]);

  end % if redraw_topo

  set(gca, 'yTick', [])

  ax(1) = min(layoutcomp.pos(:,1) - 2*layoutcomp.width);
  ax(2) = max(opt.layouttime.pos(:,1) + opt.layouttime.width/2);
  ax(3) = min(opt.layouttime.pos(:,2) - opt.layouttime.height/2);
  ax(4) = max(opt.layouttime.pos(:,2) + opt.layouttime.height/2);
  % add white space to bottom and top so channels are not out-of-axis for the majority
  % NOTE: there is another spot above with the same code, which should be kept the same as this
  % determine amount of vertical padding using cfg.verticalpadding
  if ~isnumeric(cfg.verticalpadding) && strcmp(cfg.verticalpadding, 'auto')
    % determine amount of padding using the number of channels
    if numel(cfg.channel)<=6
      wsfac = 0;
    elseif numel(cfg.channel)>6 && numel(cfg.channel)<=10
      wsfac = 0.01 *  (ax(4)-ax(3));
    else
      wsfac = 0.02 *  (ax(4)-ax(3));
    end
  else
    wsfac = cfg.verticalpadding * (ax(4)-ax(3));
  end
  ax(3) = ax(3) - wsfac;
  ax(4) = ax(4) + wsfac;
  axis(ax)
end % plotting topographies

startim = tim(1);
if nsamplepad>0
  endtim = tim(end-nsamplepad);
else
  endtim = tim(end);
end

if ~strcmp(opt.trialviewtype, 'trialsegment')
  str = sprintf('%s %d/%d, time from %g to %g s', opt.trialviewtype, opt.trlop, size(opt.trlvis,1), startim, endtim);
else
  str = sprintf('trial %d/%d: segment: %d/%d , time from %g to %g s', opt.trllock, size(opt.trlorg,1), opt.trlop, size(opt.trlvis,1), startim, endtim);
end
title(str);

% possibly adds some responsiveness if the 'thing' is clogged
drawnow

% restore the originally selected list of channels, in case a set of
% additional fixed channels was selected
cfg.channel = userchan;

setappdata(h, 'opt', opt);
setappdata(h, 'cfg', cfg);
enable_callbacks(h, true);
end % function redraw_cb
