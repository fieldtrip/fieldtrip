function [cfg] = ft_databrowser(cfg, data)

% FT_DATABROWSER can be used for visual inspection of data. Artifacts that
% were detected by artifact functions (see FT_ARTIFACT_xxx functions where
% xxx is the type of artifact) are marked. Additionally data pieces can be
% marked and unmarked as artifact by manual selection. The output cfg
% contains the updated specification of the artifacts. When visually selection
% data, a right-click will bring up a context-menu containing functions to be
% executed on the selected data. You can use your own function using cfg.selfun
% and cfg.selcfg. You can use multiple functions by giving the names/cfgs as a
% cell-array.
%
% Use as
%   cfg = ft_databrowser(cfg)
% where the configuration structure contains the reference to the dataset
% on your hard disk (see below), or use as
%   cfg = ft_databrowser(cfg, data)
% where the input data is a structure as obtained from FT_PREPROCESSING or
% from FT_COMPONENTANALYSIS.
%
% If you want to browse data that is on disk, you have to specify
%   cfg.dataset                 = string with the filename
% Instead of specifying the dataset, you can also explicitely specify the
% name of the file containing the header information and the name of the
% file containing the data, using
%   cfg.datafile                = string with the filename
%   cfg.headerfile              = string with the filename
%
% The following configuration options are supported:
%   cfg.ylim                    = vertical scaling, can be 'maxmin', 'maxabs' or [ymin ymax] (default = 'maxabs')
%   cfg.zlim                    = color scaling to apply to component topographies, 'minmax', 'maxabs' (default = 'maxmin')
%   cfg.blocksize               = duration in seconds for cutting the data up
%   cfg.trl                     = structure that defines the data segments of interest, only applicable for trial-based data
%   cfg.continuous              = 'yes' or 'no' whether the data should be interpreted as continuous or trial-based
%   cfg.channel                 = cell-array with channel labels, see FT_CHANNELSELECTION
%   cfg.plotlabels              = 'yes' (default), 'no', 'some'; whether
%                                 to plot channel labels in vertical viewmode ('some' plots one in every ten
%                                 labels; useful when plotting a large number of channels at a time)
%   cfg.ploteventlabels         = 'type=value', 'colorvalue' (default = 'type=value');
%   cfg.viewmode                = string, 'butterfly', 'vertical', 'component' for visualizing components e.g. from an ICA (default is 'butterfly')
%   cfg.artfctdef.xxx.artifact  = Nx2 matrix with artifact segments see FT_ARTIFACT_xxx functions
%   cfg.selectfeature           = string, name of feature to be selected/added (default = 'visual')
%   cfg.selectmode              = 'markartifact', 'markpeakevent', 'marktroughevent' (default = 'markartifact')
%   cfg.colorgroups             = 'sequential' 'allblack' 'labelcharx' (x = xth character in label), 'chantype' or
%                                  vector with length(data/hdr.label) defining groups (default = 'sequential')
%   cfg.channelcolormap         = COLORMAP (default = customized lines map with 15 colors)
%   cfg.selfun                  = string, name of function which is evaluated using the right-click context menu
%                                  The selected data and cfg.selcfg are passed on to this function.
%   cfg.selcfg                  = configuration options for function in cfg.selfun
%   cfg.eegscale                = number, scaling to apply to the EEG channels prior to display
%   cfg.eogscale                = number, scaling to apply to the EOG channels prior to display
%   cfg.ecgscale                = number, scaling to apply to the ECG channels prior to display
%   cfg.emgscale                = number, scaling to apply to the EMG channels prior to display
%   cfg.megscale                = number, scaling to apply to the MEG channels prior to display
%   cfg.gradscale               = number, scaling to apply to the MEG gradiometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.magscale                = number, scaling to apply to the MEG magnetometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.mychanscale             = number, scaling to apply to the channels specified in cfg.mychan
%   cfg.mychan                  = Nx1 cell-array with selection of channels
%   cfg.chanscale               = Nx1 vector with scaling factors, one per channel specified in cfg.channel
%   cfg.compscale               = string, 'global' or 'local', defines whether the colormap for the topographic scaling is
%                                  applied per topography or on all visualized components (default 'global')
%
% In case of component viewmode, a layout is required. If no layout is specified, an attempt is
% made to construct one from the sensor definition that is present in the data or specified in
% the configuration.
%   cfg.layout                  = filename of the layout, see FT_PREPARE_LAYOUT
%   cfg.elec                    = structure with electrode positions, see FT_DATATYPE_SENS
%   cfg.grad                    = structure with gradiometer definition, see FT_DATATYPE_SENS
%   cfg.elecfile                = name of file containing the electrode positions, see FT_READ_SENS
%   cfg.gradfile                = name of file containing the gradiometer definition, see FT_READ_SENS
%
% The scaling to the EEG, EOG, ECG, EMG and MEG channels is optional and
% can be used to bring the absolute numbers of the different channel types
% in the same range (e.g. fT and uV). The channel types are determined from
% the input data using FT_CHANNELSELECTION.
%
% The "artifact" field in the output cfg is a Nx2 matrix comparable to the
% "trl" matrix of FT_DEFINETRIAL. The first column of which specifying the
% beginsamples of an artifact period, the second column contains the
% endsamples of the artifactperiods.
%
% NOTE for debugging: in case the databrowser crashes, use delete(gcf) to
% kill the figure.
%
% See also FT_PREPROCESSING, FT_REJECTARTIFACT, FT_ARTIFACT_EOG,
% FT_ARTIFACT_MUSCLE, FT_ARTIFACT_JUMP, FT_ARTIFACT_MANUAL,
% FT_ARTIFACT_THRESHOLD, FT_ARTIFACT_CLIP, FT_ARTIFACT_ECG,
% FT_COMPONENTANALYSIS

% Copyright (C) 2009-2011, Robert Oostenveld, Ingrid Nieuwenhuis
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

% FIXME these should be removed or documented
% cfg.preproc
% cfg.channelcolormap
% cfg.colorgroups

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug

hasdata = (nargin>1);
hascomp = hasdata && ft_datatype(data, 'comp');

% for backward compatibility
cfg = ft_checkconfig(cfg, 'unused',     {'comps', 'inputfile', 'outputfile'});
cfg = ft_checkconfig(cfg, 'renamed',    {'zscale', 'ylim'});
cfg = ft_checkconfig(cfg, 'renamedval', {'ylim', 'auto', 'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval', {'selectmode', 'mark', 'markartifact'});

% ensure that the preproc specific options are located in the cfg.preproc substructure
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});

% set the defaults
if ~isfield(cfg, 'ylim'),            cfg.ylim = 'maxabs';                 end
if ~isfield(cfg, 'artfctdef'),       cfg.artfctdef = struct;              end
if ~isfield(cfg, 'selectfeature'),   cfg.selectfeature = 'visual';        end % string or cell-array
if ~isfield(cfg, 'selectmode'),      cfg.selectmode = 'markartifact';     end
if ~isfield(cfg, 'blocksize'),       cfg.blocksize = [];                  end % now used for both continuous and non-continuous data, defaulting done below
if ~isfield(cfg, 'preproc'),         cfg.preproc = [];                    end % see preproc for options
if ~isfield(cfg, 'selfun'),          cfg.selfun = [];                     end % default functions: 'simpleFFT','multiplotER','topoplotER','topoplotVAR','movieplotER'
if ~isfield(cfg, 'selcfg'),          cfg.selcfg = [];                     end % defaulting done below, requires layouts/etc to be processed
if ~isfield(cfg, 'colorgroups'),     cfg.colorgroups = 'sequential';      end
if ~isfield(cfg, 'channelcolormap'), cfg.channelcolormap = [0.75 0 0;0 0 1;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];   end
if ~isfield(cfg, 'eegscale'),        cfg.eegscale = [];                   end
if ~isfield(cfg, 'eogscale'),        cfg.eogscale = [];                   end
if ~isfield(cfg, 'ecgscale'),        cfg.ecgscale = [];                   end
if ~isfield(cfg, 'emgscale'),        cfg.emgscale = [];                   end
if ~isfield(cfg, 'megscale'),        cfg.megscale = [];                   end
if ~isfield(cfg, 'magscale'),        cfg.magscale = [];                   end
if ~isfield(cfg, 'gradscale'),       cfg.gradscale = [];                  end
if ~isfield(cfg, 'chanscale'),       cfg.chanscale = [];                  end
if ~isfield(cfg, 'mychanscale'),     cfg.mychanscale = [];                end
if ~isfield(cfg, 'layout'),          cfg.layout = [];                     end
if ~isfield(cfg, 'plotlabels'),      cfg.plotlabels = 'yes';              end
if ~isfield(cfg, 'event'),           cfg.event = [];                      end % this only exists for backward compatibility and should not be documented
if ~isfield(cfg, 'continuous'),      cfg.continuous = [];                 end % the default is set further down in the code, conditional on the input data
if ~isfield(cfg, 'ploteventlabels'), cfg.ploteventlabels = 'type=value';  end
cfg.zlim           = ft_getopt(cfg, 'zlim',          'maxmin');
cfg.compscale      = ft_getopt(cfg, 'compscale',     'global');


if ~isfield(cfg, 'viewmode')
  % butterfly, vertical, component
  if hascomp
    cfg.viewmode = 'component';
  else
    cfg.viewmode = 'butterfly';
  end
end

if ~isempty(cfg.chanscale)
  if ~isfield(cfg,'channel')
    warning('ignoring cfg.chanscale; this should only be used when an explicit channel selection is being made');
    cfg.chanscale = [];
  elseif numel(cfg.channel) ~= numel(cfg.chanscale)
    error('cfg.chanscale should have the same number of elements as cfg.channel');
  end
  
  % make sure chanscale is a column vector, not a row vector
  if size(cfg.chanscale,2) > size(cfg.chanscale,1)
    cfg.chanscale = cfg.chanscale';
  end
end

if ~isempty(cfg.mychanscale) && ~isfield(cfg,'mychan')
  warning('ignoring cfg.mychanscale; no channels specified in cfg.mychan');
  cfg.mychanscale = [];
end

if ~isfield(cfg, 'channel'),
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
  % read or create the layout that will be used for the topoplots
  
  if ~isempty(cfg.layout)
    tmpcfg = [];
    tmpcfg.layout = cfg.layout;
    cfg.layout = ft_prepare_layout(tmpcfg);
  else
    warning('No layout specified - will try to construct one using sensor positions');
    tmpcfg = [];
    try, tmpcfg.elec = cfg.elec; end
    try, tmpcfg.grad = cfg.grad; end
    try, tmpcfg.elecfile = cfg.elecfile; end
    try, tmpcfg.gradfile = cfg.gradfile; end
    if hasdata
      cfg.layout = ft_prepare_layout(tmpcfg, data)
    else
      cfg.layout = ft_prepare_layout(tmpcfg)
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the defaults and do some preparation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if hasdata
  % save whether data came from a timelock structure
  istimelock = strcmp(ft_datatype(data),'timelock');
  
  % check if the input data is valid for this function
  data = ft_checkdata(data, 'datatype', {'raw', 'comp'}, 'feedback', 'yes', 'hassampleinfo', 'yes');
  % fetch the header from the data structure in memory
  hdr = ft_fetch_header(data);
  
  if isfield(data, 'cfg') && ~isempty(ft_findcfg(data.cfg, 'origfs'))
    % don't use the events in case the data has been resampled
    warning('the data has been resampled, not showing the events');
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
      warning('interpreting trial-based data as continous, time-axis is no longer appropriate. t(0) now corresponds to the first sample of the first trial, and t(end) to the last sample of the last trial')
    end
  end
  
  % this is how the input data is segmented
  trlorg = zeros(numel(data.trial), 3);
  trlorg(:,[1 2]) = data.sampleinfo;
  
  % recreate offset vector (databrowser depends on this for visualisation)
  for ntrl = 1:numel(data.trial)
    trlorg(ntrl,3) = time2offset(data.time{ntrl}, data.fsample);
  end
  Ntrials = size(trlorg, 1);
  
else
  % check if the input cfg is valid for this function
  cfg = ft_checkconfig(cfg, 'dataset2files', {'yes'});
  cfg = ft_checkconfig(cfg, 'required', {'headerfile', 'datafile'});
  cfg = ft_checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
  cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});
  % read the header from file
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
  
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
    event = ft_read_event(cfg.dataset);
  end
  
  cfg.channel = ft_channelselection(cfg.channel, hdr.label);
  chansel = match_str(hdr.label, cfg.channel);
  Nchans  = length(chansel);
  
  if strcmp(cfg.continuous, 'yes')
    Ntrials = 1;
  else
    Ntrials = hdr.nTrials;
  end
  
  % FIXME in case of continuous=yes the trl should be [1 hdr.nSamples*nTrials 0]
  % and a scrollbar should be used
  
  % construct trl-matrix for data from file on disk
  trlorg = zeros(Ntrials,3);
  if strcmp(cfg.continuous, 'yes')
    trlorg(1, [1 2]) = [1 hdr.nSamples*hdr.nTrials];
  else
    for k = 1:Ntrials
      trlorg(k,[1 2]) = [1 hdr.nSamples] + [hdr.nSamples hdr.nSamples] .* (k-1);
    end
  end
end % if hasdata
if strcmp(cfg.continuous,'no') && isempty(cfg.blocksize)
  cfg.blocksize = (trlorg(1,2) - trlorg(1,1)+1) ./ hdr.Fs;
elseif strcmp(cfg.continuous,'yes') && isempty(cfg.blocksize)
  cfg.blocksize = 1;
end



% FIXME make a check for the consistency of cfg.continous, cfg.blocksize, cfg.trl and the data header

if Nchans == 0
  error('no channels to display');
end

if Ntrials == 0
  error('no trials to display');
end

if ischar(cfg.selectfeature)
  % ensure that it is a cell array
  cfg.selectfeature = {cfg.selectfeature};
end
if ~isempty(cfg.selectfeature)
  for i=1:length(cfg.selectfeature)
    if ~isfield(cfg.artfctdef, cfg.selectfeature{i})
      cfg.artfctdef.(cfg.selectfeature{i})          = [];
      cfg.artfctdef.(cfg.selectfeature{i}).artifact = zeros(0,2);
    end
  end
end

% determine the vertical scaling
if ischar(cfg.ylim)
  if hasdata
    % the first trial is used to determine the vertical scaling
    dat = data.trial{1}(chansel,:);
  else
    % one second of data is read from file to determine the vertical scaling
    dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', 1, 'endsample', round(hdr.Fs), 'chanindx', chansel, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat, 'headerformat', cfg.headerformat);
  end % if hasdata
  minval = min(dat(:));
  maxval = max(dat(:));
  switch cfg.ylim
    case 'maxabs'
      maxabs   = max(abs([minval maxval]));
      scalefac = 10^(fix(log10(maxabs)));
      maxabs   = (round(maxabs / scalefac * 100) / 100) * scalefac;
      cfg.ylim = [-maxabs maxabs];
    case 'maxmin'
      cfg.ylim = [minval maxval];
    otherwise
      error('unsupported value for cfg.ylim');
  end % switch ylim
  % zoom in a bit when viemode is vertical
  if strcmp(cfg.viewmode,'vertical')
    cfg.ylim = cfg.ylim/10;
  end
else
  if (numel(cfg.ylim) ~= 2) || ~isnumeric(cfg.ylim)
    error('cfg.ylim needs to be a 1x2 vector [ymin ymax], describing the upper and lower limits')
  end
end

% determine coloring of channels
if hasdata
  labels_all = data.label;
else
  labels_all= hdr.label;
end
if size(cfg.channelcolormap,2) ~= 3
  error('cfg.channelcolormap is not valid, size should be Nx3')
end

if isnumeric(cfg.colorgroups)
  % groups defined by user
  if length(labels_all) ~= length(cfg.colorgroups)
    error('length(cfg.colorgroups) should be length(data/hdr.label)')
  end
  R = cfg.channelcolormap(:,1);
  G = cfg.channelcolormap(:,2);
  B = cfg.channelcolormap(:,3);
  chancolors = [R(cfg.colorgroups(:)) G(cfg.colorgroups(:)) B(cfg.colorgroups(:))];
  
elseif strcmp(cfg.colorgroups, 'allblack')
  chancolors = zeros(length(labels_all),3);
  
elseif strcmp(cfg.colorgroups, 'chantype')
  type = ft_chantype(labels_all);
  [tmp1 tmp2 cfg.colorgroups] = unique(type);
  fprintf('%3d colorgroups were identified\n',length(tmp1))
  R = cfg.channelcolormap(:,1);
  G = cfg.channelcolormap(:,2);
  B = cfg.channelcolormap(:,3);
  chancolors = [R(cfg.colorgroups(:)) G(cfg.colorgroups(:)) B(cfg.colorgroups(:))];
  
elseif strcmp(cfg.colorgroups(1:9), 'labelchar')
  % groups determined by xth letter of label
  labelchar_num = str2double(cfg.colorgroups(10));
  vec_letters = num2str(zeros(length(labels_all),1));
  for iChan = 1:length(labels_all)
    vec_letters(iChan) = labels_all{iChan}(labelchar_num);
  end
  [tmp1 tmp2 cfg.colorgroups] = unique(vec_letters);
  fprintf('%3d colorgroups were identified\n',length(tmp1))
  R = cfg.channelcolormap(:,1);
  G = cfg.channelcolormap(:,2);
  B = cfg.channelcolormap(:,3);
  chancolors = [R(cfg.colorgroups(:)) G(cfg.colorgroups(:)) B(cfg.colorgroups(:))];
  
elseif strcmp(cfg.colorgroups, 'sequential')
  % no grouping
  chancolors = lines(length(labels_all));
  
else
  error('do not understand cfg.colorgroups')
end

% collect the artifacts that have been detected from cfg.artfctdef.xxx.artifact
artlabel = fieldnames(cfg.artfctdef);
sel      = zeros(size(artlabel));
artifact = cell(size(artlabel));
for i=1:length(artlabel)
  sel(i) = isfield(cfg.artfctdef.(artlabel{i}), 'artifact');
  if sel(i)
    artifact{i} = cfg.artfctdef.(artlabel{i}).artifact;
    fprintf('detected %3d %s artifacts\n', size(artifact{i}, 1), artlabel{i});
  end
end

% get the subset of the artfctdef fields that seem to contain artifacts
artifact = artifact(sel==1);
artlabel = artlabel(sel==1);

if length(artlabel) > 9
  error('only up to 9 artifacts groups supported')
end

% make artdata representing all artifacts in a "raw data" format
datendsample = max(trlorg(:,2));

artdata = [];
artdata.trial{1}       = convert_event(artifact, 'boolvec', 'endsample', datendsample); % every artifact is a "channel"
artdata.time{1}        = offset2time(0, hdr.Fs, datendsample);
artdata.label          = artlabel;
artdata.fsample        = hdr.Fs;
artdata.cfg.trl        = [1 datendsample 0];

% determine amount of unique event types (for cfg.ploteventlabels)
if ~isempty(event) && isstruct(event)
  eventtypes = unique({event.type});
else
  eventtypes = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up default functions to be available in the right-click menu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cfg.selfun - labels that are presented in rightclick menu, and is appended using ft_getuserfun(..., 'browse') later on to create a function handle
% cfg.selcfg - cfgs for functions to be executed
defselfun = [];
defselcfg = [];
% simplefft
defselfun{1} = 'simpleFFT';
tmpcfg = [];
tmpcfg.chancolors = chancolors;
defselcfg{1}  = tmpcfg;
% multiplotER
defselfun{2}  = 'multiplotER';
tmpcfg = [];
tmpcfg.layout = cfg.layout;
defselcfg{2}  = tmpcfg;
% topoplotER
defselfun{3}  = 'topoplotER';
tmpcfg = [];
tmpcfg.layout = cfg.layout;
defselcfg{3}  = tmpcfg;
% topoplotVAR
defselfun{4}  = 'topoplotVAR';
tmpcfg = [];
tmpcfg.layout = cfg.layout;
defselcfg{4}  = tmpcfg;
% movieplotER
defselfun{5}  = 'movieplotER';
tmpcfg = [];
tmpcfg.layout = cfg.layout;
tmpcfg.interactive = 'yes';
defselcfg{5}  = tmpcfg;


% add defselfuns to user-specified defselfuns
if ~iscell(cfg.selfun) && ~isempty(cfg.selfun)
  cfg.selfun = {cfg.selfun};
  cfg.selfun = [cfg.selfun defselfun];
  % do the same for the cfgs
  cfg.selcfg = {cfg.selcfg}; % assume the cfg isnt a cell
  cfg.selcfg = [cfg.selcfg defselcfg];
else
  cfg.selfun = defselfun;
  cfg.selcfg = defselcfg;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the data structures used in the GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% opt represents the global data/settings, it should contain
% - the original data, epoched or continuous
% - the artifacts represented as continuous data
% - the redraw_cb settings
% - the preproc   settings
% - the select_range_cb settings (also used in keyboard_cb)

% these elements are stored inside the figure so that the callback routines can modify them
opt = [];
if hasdata
  opt.orgdata   = data;
else
  opt.orgdata   = [];      % this means that it will look in cfg.dataset
end
if strcmp(cfg.continuous, 'yes')
  opt.trialviewtype = 'segment';
else
  opt.trialviewtype = 'trial';
end
opt.artdata     = artdata;
opt.hdr         = hdr;
opt.event       = event;
opt.trlop       = 1;          % the active trial being displayed
opt.ftsel       = find(strcmp(artlabel,cfg.selectfeature)); % current artifact/feature being selected
opt.trlorg      = trlorg;
opt.fsample     = hdr.Fs;
opt.artcolors   = [0.9686 0.7608 0.7686; 0.7529 0.7098 0.9647; 0.7373 0.9725 0.6824;0.8118 0.8118 0.8118; 0.9725 0.6745 0.4784; 0.9765 0.9176 0.5686; 0.6863 1 1; 1 0.6863 1; 0 1 0.6000];
opt.chancolors  = chancolors;
opt.cleanup     = false;      % this is needed for a corrent handling if the figure is closed (either in the corner or by "q")
opt.chanindx    = [];         % this is used to check whether the component topographies need to be redrawn
opt.eventtypes  = eventtypes;
opt.eventtypescolors = [0 0 0; 1 0 0; 0 0 1; 0 1 0; 1 0 1; 0.5 0.5 0.5; 0 1 1; 1 1 0];
opt.eventtypecolorlabels = {'black', 'red', 'blue', 'green', 'cyan', 'grey', 'light blue', 'yellow'};
opt.nanpaddata  = []; % this is used to allow horizontal scaling to be constant (when looking at last segment continous data, or when looking at segmented/zoomed-out non-continous data)
opt.trllock     = []; % this is used when zooming into trial based data

% save original layout when viewmode = component
if strcmp(cfg.viewmode,'component')
  opt.layorg    = cfg.layout;
end

% determine labelling of channels
if strcmp(cfg.plotlabels, 'yes')
  opt.plotLabelFlag = 1;
elseif strcmp(cfg.plotlabels, 'some')
  opt.plotLabelFlag = 2;
else
  opt.plotLabelFlag = 0;
end

h = figure;
setappdata(h, 'opt', opt);
setappdata(h, 'cfg', cfg);

% enable custom data cursor text
dcm = datacursormode(h);
set(dcm, 'updatefcn', @datacursortext);

% set the figure window title
funcname = mfilename();
if nargin < 2
  if isfield(cfg, 'dataset')
    dataname = cfg.dataset;
  elseif isfield(cfg, 'datafile')
    dataname = cfg.datafile;
  else
    dataname = [];
  end
else
  dataname = inputname(2);
end
set(gcf, 'Name', sprintf('%d: %s: %s', gcf, funcname, join_str(', ',dataname)));
set(gcf, 'NumberTitle', 'off');

% set zoom option to on
% zoom(h,'on')
% set(zoom(h),'actionPostCallback',@zoom_drawlabels_cb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the figure and callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(h, 'KeyPressFcn',           @keyboard_cb);
set(h, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonDownFcn'});
set(h, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonUpFcn'});
set(h, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonMotionFcn'});


% make the user interface elements for the data view
uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', opt.trialviewtype, 'userdata', 't')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'leftarrow')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'rightarrow')

uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'channel','userdata', 'c')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'uparrow')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'downarrow')

uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'horizontal', 'userdata', 'h')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'shift+leftarrow')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'shift+rightarrow')

uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'vertical', 'userdata', 'v')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'shift+downarrow')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'shift+uparrow')

% legend artifacts/features
for iArt = 1:length(artlabel)
  uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', artlabel{iArt}, 'userdata', num2str(iArt), 'position', [0.91, 0.9 - ((iArt-1)*0.09), 0.08, 0.04], 'backgroundcolor', opt.artcolors(iArt,:))
  uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', ['shift+' num2str(iArt)], 'position', [0.91, 0.855 - ((iArt-1)*0.09), 0.03, 0.04], 'backgroundcolor', opt.artcolors(iArt,:))
  uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', ['control+' num2str(iArt)], 'position', [0.96, 0.855 - ((iArt-1)*0.09), 0.03, 0.04], 'backgroundcolor', opt.artcolors(iArt,:))
end

if true % strcmp(cfg.viewmode, 'butterfly')
  % button to find label of nearest channel to datapoint
  uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'identify', 'userdata', 'i', 'position', [0.91, 0.1, 0.08, 0.05], 'backgroundcolor', [1 1 1])
end

% 'edit preproc'-button
uicontrol('tag', 'preproccfg', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string','preproc cfg','position', [0.91, 0.55 - ((iArt-1)*0.09), 0.08, 0.04],'callback',@preproc_cfg1_cb)



ft_uilayout(h, 'tag', 'labels',  'width', 0.10, 'height', 0.05);
ft_uilayout(h, 'tag', 'buttons', 'width', 0.05, 'height', 0.05);

ft_uilayout(h, 'tag', 'labels',     'style', 'pushbutton', 'callback', @keyboard_cb);
ft_uilayout(h, 'tag', 'buttons',    'style', 'pushbutton', 'callback', @keyboard_cb);
ft_uilayout(h, 'tag', 'artifactui', 'style', 'pushbutton', 'callback', @keyboard_cb);

ft_uilayout(h, 'tag', 'labels',  'retag', 'viewui');
ft_uilayout(h, 'tag', 'buttons', 'retag', 'viewui');
ft_uilayout(h, 'tag', 'viewui', 'BackgroundColor', [0.8 0.8 0.8], 'hpos', 'auto', 'vpos', 0);

definetrial_cb(h);
redraw_cb(h);


% %% Scrollbar
%
% % set initial scrollbar value
% dx = maxtime;
%
% % set scrollbar position
% fig_pos=get(gca,'position');
% scroll_pos=[fig_pos(1) fig_pos(2) fig_pos(3) 0.02];
%
% % define callback
% S=['set(gca,''xlim'',get(gcbo,''value'')+[ ' num2str(mintime) ',' num2str(maxtime) '])'];
%
% % Creating Uicontrol
% s=uicontrol('style','slider',...
%     'units','normalized','position',scroll_pos,...
%     'callback',S,'min',0,'max',0, ...
%     'visible', 'off'); %'value', xmin

% set initial scrollbar value
% dx = maxtime;
%
% % set scrollbar position
% fig_pos=get(gca,'position');
% scroll_pos=[fig_pos(1) fig_pos(2) fig_pos(3) 0.02];
%
% % define callback
% S=['set(gca,''xlim'',get(gcbo,''value'')+[ ' num2str(mintime) ',' num2str(maxtime) '])'];
%
% % Creating Uicontrol
% s=uicontrol('style','slider',...
%     'units','normalized','position',scroll_pos,...
%     'callback',S,'min',0,'max',0, ...
%     'visible', 'off'); %'value', xmin
%initialize postion of plot
% set(gca,'xlim',[xmin xmin+dx]);

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
    cfg.artfctdef.(opt.artdata.label{i}).artifact = convert_event(opt.artdata.trial{1}(i,:), 'artifact');
  end
  
  % add the updated preproc to the output
  try
    browsecfg = getappdata(h, 'cfg');
    cfg.preproc = browsecfg.preproc;
  end
  
  % add the update event to the output cfg
  cfg.event = opt.event;
  
  % do the general cleanup and bookkeeping at the end of the function
  ft_postamble debug
  ft_postamble trackconfig
  ft_postamble provenance
  ft_postamble previous data
  
end % if nargout

end % main function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanup_cb(h, eventdata)
opt = getappdata(h, 'opt');
opt.cleanup = true;
setappdata(h, 'opt', opt);
uiresume
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function definetrial_cb(h, eventdata)
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
    set(findobj(get(h,'children'),'string','trial'),'string','segment');
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
    set(findobj(get(h,'children'),'string','trialsegment'),'string',opt.trialviewtype);
    
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
    offset = begsamples - repmat(begsamples(1),[1 numel(begsamples)]); % offset for all segments compared to the first
    offset = offset + opt.trlorg(1,3);
    trlvis(:,3) = offset;
  else
    offset = begsamples - repmat(begsamples(1),[1 numel(begsamples)]);
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function help_cb(h, eventdata)
fprintf('------------------------------------------------------------------------------------\n')
fprintf('You can use the following keyboard buttons in the databrowser\n')
fprintf('1-9                : select artifact type 1-9\n');
fprintf('shift 1-9          : select previous artifact of type 1-9\n');
fprintf('                     (does not work with numpad keys)\n');
fprintf('control 1-9        : select next artifact of type 1-9\n');
fprintf('alt 1-9            : select next artifact of type 1-9\n');
fprintf('arrow-left         : previous trial\n');
fprintf('arrow-right        : next trial\n');
fprintf('shift arrow-up     : increase vertical scaling\n');
fprintf('shift arrow-down   : decrease vertical scaling\n');
fprintf('shift arrow-left   : increase horizontal scaling\n');
fprintf('shift arrow-down   : decrease horizontal scaling\n');
fprintf('s                  : toggles between cfg.selectmode options\n');
fprintf('q                  : quit\n');
fprintf('------------------------------------------------------------------------------------\n')
fprintf('\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_range_cb(h, range, cmenulab) %range 1X4 in sec relative to current trial
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

% the range should be in the displayed box
range(1) = max(opt.hpos-opt.width/2, range(1));
range(2) = max(opt.hpos-opt.width/2, range(2));
range(1) = min(opt.hpos+opt.width/2, range(1));
range(2) = min(opt.hpos+opt.width/2, range(2));
range = (range-(opt.hpos-opt.width/2)) / opt.width; % left side of the box becomes 0, right side becomes 1
range = range * (opt.hlim(2) - opt.hlim(1)) + opt.hlim(1);   % 0 becomes hlim(1), 1 becomes hlim(2)

begsample = opt.trlvis(opt.trlop,1);
endsample = opt.trlvis(opt.trlop,2);
offset    = opt.trlvis(opt.trlop,3);

% determine the selection
begsel = round(range(1)*opt.fsample+begsample-offset-1);
endsel = round(range(2)*opt.fsample+begsample-offset);
% artifact selection is now always based on begsample/endsample/offset
% -roevdmei

% the selection should always be confined to the current trial
begsel = max(begsample, begsel);
endsel = min(endsample, endsel);

% mark or execute selfun
if isempty(cmenulab)
  % the left button was clicked INSIDE a selected range, update the artifact definition or event
  
  if strcmp(cfg.selectmode, 'markartifact')
    % mark or unmark artifacts
    artval = opt.artdata.trial{1}(opt.ftsel, begsel:endsel);
    artval = any(artval,1);
    if any(artval)
      fprintf('there is overlap with the active artifact (%s), disabling this artifact\n',opt.artdata.label{opt.ftsel});
      opt.artdata.trial{1}(opt.ftsel, begsel:endsel) = 0;
    else
      fprintf('there is no overlap with the active artifact (%s), marking this as a new artifact\n',opt.artdata.label{opt.ftsel});
      opt.artdata.trial{1}(opt.ftsel, begsel:endsel) = 1;
    end
    
    % redraw only when marking (so the focus doesn't go back to the databrowser after calling selfuns
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    redraw_cb(h);
  elseif strcmp(cfg.selectmode, 'markpeakevent') || strcmp(cfg.selectmode, 'marktroughevent')
    %mark or unmark events, marking at peak/trough of window
    if any(intersect(begsel:endsel, [opt.event.sample]))
      fprintf('there is overlap with one or more event(s), disabling this/these event(s)\n');
      ind_rem = intersect(begsel:endsel, [opt.event.sample]);
      for iRemove = 1:length(ind_rem)
        opt.event([opt.event.sample]==ind_rem(iRemove)) = [];
      end
    else
      fprintf('there is no overlap with any event, adding an event to the peak/trough value\n');
      % check if only 1 chan, other wise not clear max in which channel. %
      % ingnie: would be cool to add the option to select the channel when multiple channels
      if size(opt.curdat.trial{1},1) > 1
        error('cfg.selectmode = ''markpeakevent'' and ''marktroughevent'' only supported with 1 channel in the data')
      end
      if strcmp(cfg.selectmode, 'markpeakevent')
        [dum ind_minmax] = max(opt.curdat.trial{1}(begsel-begsample+1:endsel-begsample+1));
        val = 'peak';
      elseif strcmp(cfg.selectmode, 'marktroughevent')
        [dum ind_minmax] = min(opt.curdat.trial{1}(begsel-begsample+1:endsel-begsample+1));
        val = 'trough';
      end
      samp_minmax = begsel + ind_minmax - 1;
      event_new.type     = 'ft_databrowser_manual';
      event_new.sample   = samp_minmax;
      event_new.value    = val;
      event_new.duration = 1;
      event_new.offset   = 0;
      % add new event to end opt.event
      % check if events are in order now
      if  min(diff([opt.event.sample]))>0
        % add new event in line with old ones
        nearest_event = nearest([opt.event.sample], samp_minmax);
        if opt.event(nearest_event).sample > samp_minmax
          %place new event before nearest
          ind_event_new = nearest_event;
        else
          %place new event after nearest
          ind_event_new = nearest_event +1;
        end
        event_lastpart = opt.event(ind_event_new:end);
        opt.event(ind_event_new) = event_new;
        opt.event(ind_event_new+1:end+1) = event_lastpart;
      else
        %just add to end
        opt.event(end+1) = event_new;
      end
      clear event_new ind_event_new event_lastpart val dum ind_minmax
    end
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
  seldata.label    = opt.curdat.label;
  seldata.time{1}  = offset2time(offset+begsel-begsample, opt.fsample, endsel-begsel+1);
  seldata.trial{1} = ft_fetch_data(opt.curdat, 'begsample', begsel, 'endsample', endsel);
  seldata.fsample  = opt.fsample;
  seldata.cfg.trl  = [begsel endsel offset];
  
  % prepare input
  funhandle = ft_getuserfun(cmenulab,'browse');
  funcfg = cfg.selcfg{selfunind};
  % get windowname and give as input (can be used for the other functions as well, not implemented yet)
  if ~strcmp(opt.trialviewtype,'trialsegment')
    str = sprintf('%s %d/%d, time from %g to %g s', opt.trialviewtype, opt.trlop, size(opt.trlvis,1), seldata.time{1}(1), seldata.time{1}(end));
  else
    str = sprintf('trial %d/%d: segment: %d/%d , time from %g to %g s', opt.trllock, size(opt.trlorg,1), opt.trlop, size(opt.trlvis,1), seldata.time{1}(1), seldata.time{1}(end));
  end
  funcfg.figurename = [cmenulab ': ' str];
  feval(funhandle, funcfg, seldata);
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function preproc_cfg1_cb(h,eventdata)
parent = get(h,'parent');
cfg = getappdata(parent, 'cfg');

% parse cfg.preproc
if ~isempty(cfg.preproc)
  tmpcfg = cfg.preproc;
  cfg = [];
  cfg.preproc = tmpcfg;
  code = printstruct('cfg', cfg);
else
  code = [];
end

% add descriptive lines
sep     = sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
descrip = sprintf('%% Add/change cfg options for on-the-fly preprocessing\n%% Use as cfg.preproc.xxx\n');
code    = [sep descrip sep code];


% make figure displaying the edit box
pph = figure;
axis off
% add save button
uicontrol('tag', 'preproccfg_l2', 'parent', pph, 'units', 'normalized', 'style', 'pushbutton', 'string','save and close','position', [0.81, 0.6 , 0.18, 0.10],'callback',@preproc_cfg2_cb);

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
set(ppeh, 'FontName', 'Courier');
set(ppeh, 'FontSize', 12);
set(ppeh, 'string', code);


% add handle for the edit style to figure
setappdata(pph,'superparent', parent); % superparent is the main ft_databrowser window
setappdata(pph,'ppeh', ppeh);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function preproc_cfg2_cb(h,eventdata)
parent = get(h,'parent');
superparent = getappdata(parent,'superparent');
ppeh   = getappdata(parent,'ppeh');
code = get(ppeh, 'string');

% remove descriptive lines (so they don't display on command line)
code = cellstr(code(5:end,:));
% get rid of empty lines and white space
remind = [];
for iline = 1:numel(code)
  code{iline} = strtrim(code{iline});
  if isempty(code{iline})
    remind = [remind iline];
  end
end
code(remind) = [];

if ~isempty(code)
  ispreproccfg = strncmp(code,'cfg.preproc.',12);
  if ~all(ispreproccfg)
    errordlg('cfg-options must be specified as cfg.preproc.xxx','cfg.preproc editor','modal')
  end
  % eval the code
  for icomm = 1:numel(code)
    eval([code{icomm} ';']);
  end
  
  % check for cfg and output into the original appdata-window
  if ~exist('cfg','var')
    cfg = [];
    cfg.preproc = [];
  end
  maincfg = getappdata(superparent,'cfg');
  maincfg.preproc = cfg.preproc;
  setappdata(superparent,'cfg',maincfg)
end

close(parent)
redraw_cb(superparent)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keyboard_cb(h, eventdata)

if isempty(eventdata)
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parseKeyboardEvent(eventdata);
end
% get focus back to figure
if ~strcmp(get(h, 'type'), 'figure')
  set(h, 'enable', 'off');
  drawnow;
  set(h, 'enable', 'on');
end

h = getparent(h);
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

switch key
  case {'1' '2' '3' '4' '5' '6' '7' '8' '9'}
    % switch to another artifact type
    opt.ftsel = str2double(key);
    numart = size(opt.artdata.trial{1}, 1);
    if opt.ftsel > numart
      fprintf('data has no artifact type %i \n', opt.ftsel)
    else
      setappdata(h, 'opt', opt);
      setappdata(h, 'cfg', cfg);
      fprintf('switching to the "%s" artifact\n', opt.artdata.label{opt.ftsel});
      redraw_cb(h, eventdata);
    end
  case {'shift+1' 'shift+2' 'shift+3' 'shift+4' 'shift+5' 'shift+6' 'shift+7' 'shift+8' 'shift+9'}
    % go to previous artifact
    opt.ftsel = str2double(key(end));
    numart = size(opt.artdata.trial{1}, 1);
    if opt.ftsel > numart
      fprintf('data has no artifact type %i \n', opt.ftsel)
    else
      cursam = opt.trlvis(opt.trlop,1);
      artsam = find(opt.artdata.trial{1}(opt.ftsel,1:cursam-1), 1, 'last');
      if isempty(artsam)
        fprintf('no earlier "%s" artifact found\n', opt.artdata.label{opt.ftsel});
      else
        fprintf('going to previous "%s" artifact\n', opt.artdata.label{opt.ftsel});
        if opt.trlvis(nearest(opt.trlvis(:,1),artsam),1) < artsam
          arttrl = nearest(opt.trlvis(:,1),artsam);
        else
          arttrl = nearest(opt.trlvis(:,1),artsam)-1;
        end
        opt.trlop = arttrl;
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h, eventdata);
      end
    end
  case {'control+1' 'control+2' 'control+3' 'control+4' 'control+5' 'control+6' 'control+7' 'control+8' 'control+9' 'alt+1' 'alt+2' 'alt+3' 'alt+4' 'alt+5' 'alt+6' 'alt+7' 'alt+8' 'alt+9'}
    % go to next artifact
    opt.ftsel = str2double(key(end));
    numart = size(opt.artdata.trial{1}, 1);
    if opt.ftsel > numart
      fprintf('data has no artifact type %i \n', opt.ftsel)
    else
      cursam = opt.trlvis(opt.trlop,2);
      artsam = find(opt.artdata.trial{1}(opt.ftsel,cursam+1:end), 1, 'first') + cursam;
      if isempty(artsam)
        fprintf('no later "%s" artifact found\n', opt.artdata.label{opt.ftsel});
      else
        fprintf('going to next "%s" artifact\n', opt.artdata.label{opt.ftsel});
        if opt.trlvis(nearest(opt.trlvis(:,1),artsam),1) < artsam
          arttrl = nearest(opt.trlvis(:,1),artsam);
        else
          arttrl = nearest(opt.trlvis(:,1),artsam)-1;
        end
        opt.trlop = arttrl;
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
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    delete(findobj(h,'tag', 'chanlabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
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
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    delete(findobj(h,'tag', 'chanlabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
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
  case 'q'
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    cleanup_cb(h);
  case 't'
    % select the trial to display
    if ~strcmp(opt.trialviewtype,'trialsegment')
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
      response = ['[' response{1} ']']; % convert to string and add brackets, just to ensure that str2num will work
      if strcmp(response, '[maxmin]')
        minval = min(opt.curdat.trial{1}(:));
        maxval = max(opt.curdat.trial{1}(:));
        cfg.ylim = [minval maxval];
      elseif strcmp(response, '[maxabs]')
        minval = min(opt.curdat.trial{1}(:));
        maxval = max(opt.curdat.trial{1}(:));
        cfg.ylim = [-max(abs([minval maxval])) max(abs([minval maxval]))];
      else
        tmp = str2num(response);
        if numel(tmp)==2
          cfg.ylim = tmp;
        else
          warning('incorrect specification of cfg.ylim, not changing the limits for the vertical axes')
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
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    delete(findobj(h,'tag', 'chanlabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
    redraw_cb(h, eventdata);
  case 'i'
    if strcmp(cfg.viewmode, 'butterfly')
      delete(findobj(h,'tag', 'identify'));
      % click in data and get name of nearest channel
      fprintf('click in the figure to identify the name of the closest channel\n');
      val = ginput(1);
      pos = val(1);
      % transform 'val' to match data
      val(1) = val(1) * range(opt.hlim) + opt.hlim(1);
      val(2) = val(2) * range(opt.vlim) + opt.vlim(1);
      channame = val2nearestchan(opt.curdat,val);
      channb = match_str(opt.curdat.label,channame);
      fprintf('channel name: %s\n',channame);
      redraw_cb(h, eventdata);
      ft_plot_text(pos, 0.9, channame, 'FontSize', 16, 'tag', 'identify','interpreter','none');
      if ~ishold
        hold on
        ft_plot_vector(opt.curdat.time{1}, opt.curdat.trial{1}(channb,:), 'box', false, 'tag', 'identify', ...
          'hpos', opt.laytime.pos(1,1), 'vpos', opt.laytime.pos(1,2), 'width', opt.laytime.width(1), 'height', opt.laytime.height(1), 'hlim', opt.hlim, 'vlim', opt.vlim, ...
          'color', 'k', 'linewidth', 2);
        hold off
      else
        ft_plot_vector(opt.curdat.time{1}, opt.curdat.trial{1}(channb,:), 'box', false, 'tag', 'identify', ...
          'hpos', opt.laytime.pos(1,1), 'vpos', opt.laytime.pos(1,2), 'width', opt.laytime.width(1), 'height', opt.laytime.height(1), 'hlim', opt.hlim, 'vlim', opt.vlim, ...
          'color', 'k', 'linewidth', 2);
      end
    else
      warning('only supported with cfg.viewmode=''butterfly''');
    end
  case 's'
    % toggle between selectmode options: switch from 'markartifact', to 'markpeakevent' to 'marktroughevent' and back with on screen feedback
    curstate = find(strcmp(cfg.selectmode, {'markartifact', 'markpeakevent', 'marktroughevent'}));
    if curstate == 1
      cfg.selectmode = 'markpeakevent';
    elseif curstate == 2
      cfg.selectmode = 'marktroughevent';
    elseif curstate == 3
      cfg.selectmode = 'markartifact';
    end
    fprintf('switching to selectmode = %s\n',cfg.selectmode);
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    redraw_cb(h, eventdata);
  case 'control+control'
    % do nothing
  case 'shift+shift'
    % do nothing
  case 'alt+alt'
    % do nothing
  otherwise
    setappdata(h, 'opt', opt);
    setappdata(h, 'cfg', cfg);
    help_cb(h);
end
uiresume(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toggle_viewmode_cb(h, eventdata, varargin)
% FIXME should be used
opt = guidata(getparent(h));
if ~isempty(varargin) && ischar(varargin{1})
  cfg.viewmode = varargin{1};
end
guidata(getparent(h), opt);
delete(findobj(h,'tag', 'chanlabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
redraw_cb(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function redraw_cb(h, eventdata)
h = getparent(h);
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

figure(h); % ensure that the calling figure is in the front

%fprintf('redrawing with viewmode %s\n', cfg.viewmode);

begsample = opt.trlvis(opt.trlop, 1);
endsample = opt.trlvis(opt.trlop, 2);
offset    = opt.trlvis(opt.trlop, 3);
chanindx  = match_str(opt.hdr.label, cfg.channel);

if ~isempty(opt.event) && isstruct(opt.event)
  % select only the events in the current time window
  event     = opt.event;
  evtsample = [event(:).sample];
  event     = event(evtsample>=begsample & evtsample<=endsample);
else
  event = [];
end

if isempty(opt.orgdata)
  dat = ft_read_data(cfg.datafile, 'header', opt.hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat, 'headerformat', cfg.headerformat);
else
  dat = ft_fetch_data(opt.orgdata, 'header', opt.hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'allowoverlap', true); % ALLOWING OVERLAPPING TRIALS
end
art = ft_fetch_data(opt.artdata, 'begsample', begsample, 'endsample', endsample);

% apply preprocessing and determine the time axis
[dat, lab, tim] = preproc(dat, opt.hdr.label(chanindx), offset2time(offset, opt.fsample, size(dat,2)), cfg.preproc);

% add NaNs to data for plotting purposes. NaNs will be added when requested horizontal scaling is longer than the data.
nsamplepad = opt.nanpaddata(opt.trlop);
if nsamplepad>0
  dat = [dat NaN(numel(lab), opt.nanpaddata(opt.trlop))];
  tim = [tim linspace(tim(end),tim(end)+nsamplepad*mean(diff(tim)),nsamplepad)];  % possible machine precision error here
end
opt.curdat.label      = lab;
opt.curdat.time{1}    = tim;
opt.curdat.trial{1}   = dat;
opt.curdat.fsample    = opt.fsample;
opt.curdat.sampleinfo = [begsample endsample offset];

% apply scaling to selected channels
% using wildcard to support subselection of channels
if ~isempty(cfg.eegscale)
  chansel = match_str(lab, ft_channelselection('EEG', lab));
  dat(chansel,:) = dat(chansel,:) .* cfg.eegscale;
end
if ~isempty(cfg.eogscale)
  chansel = match_str(lab, ft_channelselection('EOG', lab));
  dat(chansel,:) = dat(chansel,:) .* cfg.eogscale;
end
if ~isempty(cfg.ecgscale)
  chansel = match_str(lab, ft_channelselection('ECG', lab));
  dat(chansel,:) = dat(chansel,:) .* cfg.ecgscale;
end
if ~isempty(cfg.emgscale)
  chansel = match_str(lab, ft_channelselection('EMG', lab));
  dat(chansel,:) = dat(chansel,:) .* cfg.emgscale;
end
if ~isempty(cfg.megscale)
  type = opt.hdr.grad.type;
  chansel = match_str(lab, ft_channelselection('MEG', lab, type));
  dat(chansel,:) = dat(chansel,:) .* cfg.megscale;
end
if ~isempty(cfg.magscale)
  chansel = match_str(lab, ft_channelselection('MEGMAG', lab));
  dat(chansel,:) = dat(chansel,:) .* cfg.magscale;
end
if ~isempty(cfg.gradscale)
  chansel = match_str(lab, ft_channelselection('MEGGRAD', lab));
  dat(chansel,:) = dat(chansel,:) .* cfg.gradscale;
end
if ~isempty(cfg.chanscale)
  chansel = match_str(lab, ft_channelselection(cfg.channel, lab));
  dat(chansel,:) = dat(chansel,:) .* repmat(cfg.chanscale,1,size(dat,2));
end
if ~isempty(cfg.mychanscale)
  chansel = match_str(lab, ft_channelselection(cfg.mychan, lab));
  dat(chansel,:) = dat(chansel,:) .* cfg.mychanscale;
end

% to assure current feature is plotted on top
ordervec = 1:length(opt.artdata.label);
ordervec(opt.ftsel) = [];
ordervec(end+1) = opt.ftsel;

% FIXME speedup ft_prepare_layout
if strcmp(cfg.viewmode, 'butterfly')
  laytime = [];
  laytime.label = {'dummy'};
  laytime.pos = [0.5 0.5];
  laytime.width = 1;
  laytime.height = 1;
  opt.laytime = laytime;
else
  % this needs to be reconstructed if the channel selection changes
  tmpcfg = [];
  if strcmp(cfg.viewmode, 'component')
    tmpcfg.layout  = 'vertical';
  else
    tmpcfg.layout  = cfg.viewmode;
  end
  tmpcfg.channel = cfg.channel;
  tmpcfg.skipcomnt = 'yes';
  tmpcfg.skipscale = 'yes';
  tmpcfg.showcallinfo = 'no';
  opt.laytime = ft_prepare_layout(tmpcfg, opt.orgdata);
end

% determine the position of the channel/component labels relative to the real axes
% FIXME needs a shift to the left for components
labelx = opt.laytime.pos(:,1) - opt.laytime.width/2 - 0.01;
labely = opt.laytime.pos(:,2);

% determine the total extent of all virtual axes relative to the real axes
ax(1) = min(opt.laytime.pos(:,1) - opt.laytime.width/2);
ax(2) = max(opt.laytime.pos(:,1) + opt.laytime.width/2);
ax(3) = min(opt.laytime.pos(:,2) - opt.laytime.height/2);
ax(4) = max(opt.laytime.pos(:,2) + opt.laytime.height/2);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('plotting artifacts...\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete(findobj(h, 'tag', 'artifact'));

for j = ordervec
  tmp = diff([0 art(j,:) 0]);
  artbeg = find(tmp==+1);
  artend = find(tmp==-1) - 1;
  
  for k=1:numel(artbeg)
    xpos = [tim(artbeg(k)) tim(artend(k))] + ([-.5 +.5]./opt.fsample);
    h_artifact = ft_plot_box([xpos -1 1], 'facecolor', opt.artcolors(j,:), 'edgecolor', 'none', 'tag', 'artifact',  ...
      'hpos', opt.hpos, 'vpos', opt.vpos, 'width', opt.width, 'height', opt.height, 'hlim', opt.hlim, 'vlim', [-1 1]);
  end
end % for each of the artifact channels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('plotting events...\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete(findobj(h, 'tag', 'event'));

if any(strcmp(cfg.viewmode, {'butterfly', 'component', 'vertical'}))
  
  if strcmp(cfg.ploteventlabels , 'colorvalue') && ~isempty(opt.event)
    eventlabellegend = [];
    for iType = 1:length(opt.eventtypes)
      eventlabellegend = [eventlabellegend sprintf('%s = %s\n',opt.eventtypes{iType},opt.eventtypecolorlabels{iType})];
    end
    fprintf(eventlabellegend);
  end
  
  % save stuff to able to shift event labels downwards when they occur at the same time-point
  eventcol = cell(1,numel(event));
  eventstr = cell(1,numel(event));
  eventtim = NaN(1,numel(event));
  
  % gather event info and plot lines
  for ievent = 1:numel(event)
    try
      if strcmp(cfg.ploteventlabels , 'type=value')
        if isempty(event(ievent).value)
          eventstr{ievent} = '';
        else
          eventstr{ievent} = sprintf('%s = %s', event(ievent).type, num2str(event(ievent).value)); % value can be both number and string
        end
        eventcol{ievent} = 'k';
      elseif strcmp(cfg.ploteventlabels , 'colorvalue')
        eventcol{ievent} = opt.eventtypescolors(match_str(opt.eventtypes, event(ievent).type),:);
        eventstr{ievent} = sprintf('%s', num2str(event(ievent).value)); % value can be both number and string
      end
    catch
      eventstr{ievent} = 'unknown';
      eventcol{ievent} = 'k';
    end
    eventtim(ievent) = (event(ievent).sample-begsample)/opt.fsample + opt.hlim(1);
    
    lh = ft_plot_line([eventtim(ievent) eventtim(ievent)], [-1 1], 'tag', 'event', 'color', eventcol{ievent}, ...
      'hpos', opt.hpos, 'vpos', opt.vpos, 'width', opt.width, 'height', opt.height, 'hlim', opt.hlim, 'vlim', [-1 1]);
    
    % store this data in the line object so that it can be displayed in the
    % data cursor (see subfunction datacursortext below)
    setappdata(lh, 'ft_databrowser_linetype', 'event');
    setappdata(lh, 'ft_databrowser_eventtime', eventtim(ievent));
    setappdata(lh, 'ft_databrowser_eventtype', event(ievent).type);
    setappdata(lh, 'ft_databrowser_eventvalue', event(ievent).value);
  end
  % count the consecutive occurrence of each time point
  concount = NaN(1,numel(event));
  for ievent = 1:numel(event)
    concount(ievent) = sum(eventtim(ievent)==eventtim(1:ievent-1));
  end
  % plot event labels
  for ievent = 1:numel(event)
    ft_plot_text(eventtim(ievent), 0.9-concount(ievent)*.06, eventstr{ievent}, 'tag', 'event', 'Color', eventcol{ievent}, ...
      'hpos', opt.hpos, 'vpos', opt.vpos, 'width', opt.width, 'height', opt.height, 'hlim', opt.hlim, 'vlim', [-1 1],'interpreter','none');
  end
  
end % if viewmode appropriate for events

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('plotting data...\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete(findobj(h,'tag', 'timecourse'));
delete(findobj(h,'tag', 'identify'));

% not removing channel labels, they cause the bulk of redrawing time for the slow text function (note, interpreter = none hardly helps)
% warning, when deleting the labels using the line below, one can easily tripple the excution time of redrawing in viewmode = vertical (see bug 2065)
%delete(findobj(h, 'tag', 'chanlabel'));

if strcmp(cfg.viewmode, 'butterfly')
  set(gca,'ColorOrder',opt.chancolors(chanindx,:)) % plot vector does not clear axis, therefore this is possible
  ft_plot_vector(tim, dat, 'box', false, 'tag', 'timecourse', ...
    'hpos', opt.laytime.pos(1,1), 'vpos', opt.laytime.pos(1,2), 'width', opt.laytime.width(1), 'height', opt.laytime.height(1), 'hlim', opt.hlim, 'vlim', opt.vlim);
  
  
  % two ticks per channel
  yTick = sort([opt.laytime.pos(:,2)+(opt.laytime.height/2); ...
    opt.laytime.pos(:,2)+(opt.laytime.height/4); ...
    opt.laytime.pos(:,2);                        ...
    opt.laytime.pos(:,2)-(opt.laytime.height/4); ...
    opt.laytime.pos(:,2)-(opt.laytime.height/2)]);
  
  yTickLabel = {num2str(yTick.*range(opt.vlim) + opt.vlim(1))};
  
  set(gca, 'yTick', yTick);
  set(gca, 'yTickLabel', yTickLabel)
  
elseif any(strcmp(cfg.viewmode, {'component', 'vertical'}))
  
    % determine channel indices into data outside of loop
  laysels = match_str(opt.laytime.label, opt.hdr.label);
  
  for i = 1:length(chanindx)
    if strcmp(cfg.viewmode, 'component')
      color = 'k';
    else
      color = opt.chancolors(chanindx(i),:);
    end
    datsel = i;
    laysel = laysels(i);

    if ~isempty(datsel) && ~isempty(laysel)
      % only plot labels when current chanlabel objects are less then the total number of channels (see bug 2065)
      % this is a cheap quick fix. If it causes error in plotting components, do this conditional on viewmode
      if numel(findobj(h,'tag', 'chanlabel'))<numel(chanindx)
        if opt.plotLabelFlag == 1 || (opt.plotLabelFlag == 2 && mod(i,10)==0)
          ft_plot_text(labelx(laysel), labely(laysel), opt.hdr.label(chanindx(i)), 'tag', 'chanlabel', 'HorizontalAlignment', 'right','interpreter','none','FontUnits','normalized','FontSize',8);
        end
      end
      
      lh = ft_plot_vector(tim, dat(datsel, :), 'box', false, 'color', color, 'tag', 'timecourse', ...
        'hpos', opt.laytime.pos(laysel,1), 'vpos', opt.laytime.pos(laysel,2), 'width', opt.laytime.width(laysel), 'height', opt.laytime.height(laysel), 'hlim', opt.hlim, 'vlim', opt.vlim);
      
      % store this data in the line object so that it can be displayed in the
      % data cursor (see subfunction datacursortext below)
      setappdata(lh, 'ft_databrowser_linetype', 'channel');
      setappdata(lh, 'ft_databrowser_label', opt.hdr.label(chanindx(i)));
      setappdata(lh, 'ft_databrowser_xaxis', tim);
      setappdata(lh, 'ft_databrowser_yaxis', dat(datsel,:));
    end
  end
  
  if length(chanindx)>19
    % no space for yticks
    yTick = [];
    yTickLabel = [];
  elseif length(chanindx)> 6
    % one tick per channel
    yTick = sort([opt.laytime.pos(:,2)+(opt.laytime.height(laysel)/4); ...
      opt.laytime.pos(:,2)-(opt.laytime.height(laysel)/4)]);
    yTickLabel = {[.25 .75] .* range(opt.vlim) + opt.vlim(1)};
  else
    % two ticks per channel
    yTick = sort([opt.laytime.pos(:,2)+(opt.laytime.height(laysel)/2); ...
      opt.laytime.pos(:,2)+(opt.laytime.height(laysel)/4); ...
      opt.laytime.pos(:,2)-(opt.laytime.height(laysel)/4); ...
      opt.laytime.pos(:,2)-(opt.laytime.height(laysel)/2)]);
    yTickLabel = {[.0 .25 .75 1] .* range(opt.vlim) + opt.vlim(1)};
  end
  
  yTickLabel = repmat(yTickLabel, 1, length(chanindx));
  set(gca, 'yTick', yTick);
  set(gca, 'yTickLabel', yTickLabel);
  
else
  % the following is implemented for 2column, 3column, etcetera.
  % it also works for topographic layouts, such as CTF151
  
  % determine channel indices into data outside of loop
  laysels = match_str(opt.laytime.label, opt.hdr.label);
  
  for i = 1:length(chanindx)
    color = opt.chancolors(chanindx(i),:);
    datsel = i;
    laysel = laysels(i);
    
    if ~isempty(datsel) && ~isempty(laysel)
      
      lh = ft_plot_vector(tim, dat(datsel, :), 'box', false, 'color', color, 'tag', 'timecourse', ...
        'hpos', opt.laytime.pos(laysel,1), 'vpos', opt.laytime.pos(laysel,2), 'width', opt.laytime.width(laysel), 'height', opt.laytime.height(laysel), 'hlim', opt.hlim, 'vlim', opt.vlim);
      
      % store this data in the line object so that it can be displayed in the
      % data cursor (see subfunction datacursortext below)
      setappdata(lh, 'ft_databrowser_linetype', 'channel');
      setappdata(lh, 'ft_databrowser_label', opt.hdr.label(chanindx(i)));
      setappdata(lh, 'ft_databrowser_xaxis', tim);
      setappdata(lh, 'ft_databrowser_yaxis', dat(datsel,:));
    end
  end
  
  % ticks are not supported with such a layout
  yTick = [];
  yTickLabel = [];
  
  yTickLabel = repmat(yTickLabel, 1, length(chanindx));
  set(gca, 'yTick', yTick);
  set(gca, 'yTickLabel', yTickLabel);
  
end % if strcmp viewmode

if any(strcmp(cfg.viewmode, {'butterfly', 'component', 'vertical'}))
  nticks = 11;
  xTickLabel = cellstr(num2str( linspace(tim(1), tim(end), nticks)' , '%1.2f'))';
  if nsamplepad>0
    nlabindat = sum(linspace(tim(1), tim(end), nticks) < tim(end-nsamplepad));
    xTickLabel(nlabindat+1:end) = repmat({' '},[1 nticks-nlabindat]);
  end
  set(gca, 'xTick', linspace(ax(1), ax(2), nticks))
  set(gca, 'xTickLabel', xTickLabel)
  xlabel('time');
else
  set(gca, 'xTick', [])
  set(gca, 'xTickLabel', [])
end

if strcmp(cfg.viewmode, 'component')
  
  % determine the position of each of the original channels for the topgraphy
  laychan = opt.layorg;
  
  % determine the position of each of the topographies
  laytopo.pos(:,1)  = opt.laytime.pos(:,1) - opt.laytime.width/2 - opt.laytime.height;
  laytopo.pos(:,2)  = opt.laytime.pos(:,2) + opt.laytime.height/2;
  laytopo.width     = opt.laytime.height;
  laytopo.height    = opt.laytime.height;
  laytopo.label     = opt.laytime.label;
  
  if ~isequal(opt.chanindx, chanindx)
    opt.chanindx = chanindx;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fprintf('plotting component topographies...\n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delete(findobj(h,'tag', 'topography'));
    
    [sel1, sel2] = match_str(opt.orgdata.topolabel, laychan.label);
    chanx = laychan.pos(sel2,1);
    chany = laychan.pos(sel2,2);
    
    if strcmp(cfg.compscale, 'global')
      for i=1:length(chanindx) % loop through all components to get max and min
        zmin(i) = min(opt.orgdata.topo(sel1,chanindx(i)));
        zmax(i) = max(opt.orgdata.topo(sel1,chanindx(i)));
      end
      
      if strcmp(cfg.zlim, 'maxmin')
        zmin = min(zmin);
        zmax = max(zmax);
      elseif strcmp(cfg.zlim, 'maxabs')
        zmax = max([abs(zmin) abs(zmax)]);
        zmin = -zmax;
      else
        error('configuration option for component scaling could not be recognized');
      end
    end
    
    for i=1:length(chanindx)
      % plot the topography of this component
      laysel = match_str(opt.laytime.label, opt.hdr.label(chanindx(i)));
      chanz = opt.orgdata.topo(sel1,chanindx(i));
      
      if strcmp(cfg.compscale, 'local')
        % compute scaling factors here
        if strcmp(cfg.zlim, 'maxmin')
          zmin = min(chanz);
          zmax = max(chanz);
        elseif strcmp(cfg.zlim, 'maxabs')
          zmax = max(abs(chanz));
          zmin = -zmax;
        end
      end
      
      % scaling
      chanz = (chanz - zmin) ./  (zmax- zmin);
      
      % laychan is the actual topo layout, in pixel units for .mat files
      % laytopo is a vertical layout determining where to plot each topo,
      %   with one entry per component
      
      
      ft_plot_topo(chanx, chany, chanz, 'mask', ...
        laychan.mask, 'interplim', 'mask', 'outline', ...
        laychan.outline, 'tag', 'topography', ...
        'hpos', laytopo.pos(laysel,1)-laytopo.width(laysel)/2,...
        'vpos', laytopo.pos(laysel,2)-laytopo.height(laysel)/2,...
        'width', laytopo.width(laysel), 'height', laytopo.height(laysel), 'gridscale', 45);
      
      %axis equal
      %drawnow
    end
    
    caxis([0 1]);
    
  end % if redraw_topo
  
  set(gca, 'yTick', [])
  
  ax(1) = min(laytopo.pos(:,1) - laytopo.width);
  ax(2) = max(opt.laytime.pos(:,1) + opt.laytime.width/2);
  ax(3) = min(opt.laytime.pos(:,2) - opt.laytime.height/2);
  ax(4) = max(opt.laytime.pos(:,2) + opt.laytime.height/2);
  axis(ax)
end % plotting topographies

startim = tim(1);
if nsamplepad>0
  endtim = tim(end-nsamplepad);
else
  endtim = tim(end);
end

if ~strcmp(opt.trialviewtype,'trialsegment')
  str = sprintf('%s %d/%d, time from %g to %g s', opt.trialviewtype, opt.trlop, size(opt.trlvis,1), startim, endtim);
else
  str = sprintf('trial %d/%d: segment: %d/%d , time from %g to %g s', opt.trllock, size(opt.trlorg,1), opt.trlop, size(opt.trlvis,1), startim, endtim);
end
title(str);

% possibly adds some responsiveness if the 'thing' is clogged
drawnow

setappdata(h, 'opt', opt);
setappdata(h, 'cfg', cfg);
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = parseKeyboardEvent(eventdata)

key = eventdata.Key;

% handle possible numpad events (different for Windows and UNIX systems)
% NOTE: shift+numpad number does not work on UNIX, since the shift
% modifier is always sent for numpad events
if isunix()
  shiftInd = match_str(eventdata.Modifier, 'shift');
  if ~isnan(str2double(eventdata.Character)) && ~isempty(shiftInd)
    % now we now it was a numpad keystroke (numeric character sent AND
    % shift modifier present)
    key = eventdata.Character;
    eventdata.Modifier(shiftInd) = []; % strip the shift modifier
  end
elseif ispc()
  if strfind(eventdata.Key, 'numpad')
    key = eventdata.Character;
  end
end

if ~isempty(eventdata.Modifier)
  key = [eventdata.Modifier{1} '+' key];
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cursortext = datacursortext(obj, event_obj)
pos = get(event_obj, 'position');

linetype = getappdata(event_obj.Target, 'ft_databrowser_linetype');

if strcmp(linetype, 'event')
  cursortext = sprintf('%s = %d\nt = %g s',...
    getappdata(event_obj.Target, 'ft_databrowser_eventtype'),...
    getappdata(event_obj.Target, 'ft_databrowser_eventvalue'),...
    getappdata(event_obj.Target, 'ft_databrowser_eventtime'));
elseif strcmp(linetype, 'channel')
  % get plotted x axis
  plottedX = get(event_obj.Target, 'xdata');
  
  % determine values of data at real x axis
  timeAxis = getappdata(event_obj.Target, 'ft_databrowser_xaxis');
  dataAxis = getappdata(event_obj.Target, 'ft_databrowser_yaxis');
  tInd = nearest(plottedX, pos(1));
  
  % get label
  chanLabel = getappdata(event_obj.Target, 'ft_databrowser_label');
  chanLabel = chanLabel{1};
  
  cursortext = sprintf('t = %g\n%s = %g', timeAxis(tInd), chanLabel,...
    dataAxis(tInd));
else
  cursortext = '<no cursor available>';
  % explicitly tell the user there is no info because the x-axis and
  % y-axis do not correspond to real data values (both are between 0 and
  % 1 always)
end
end % function
