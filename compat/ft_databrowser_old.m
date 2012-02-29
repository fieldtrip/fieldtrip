function [cfg] = ft_databrowser_old(cfg, data)

% FT_DATABROWSER can be used for visual inspection of data. Artifacts that were detected
% by artifact functions (see FT_ARTIFACT_xxx functions where xxx is the type of artifact)
% are marked. Additionally data pieces can be marked and unmarked as artifact by
% manual selection. The output cfg contains the updated artifactdef field.
%
% Use as
%   cfg = ft_databrowser(cfg)
%   required configuration options:
%   cfg.dataset or both cfg.headerfile and cfg.datafile
% or as
%   cfg = ft_databrowser(cfg, data)
% with the data as obtained from FT_PREPROCESSING
%
% The following configuration options are supported:
%   cfg.trl                     = structure that defines the data segments of interest. See FT_DEFINETRIAL
%   cfg.continuous              = 'yes' or 'no' wh ether the file contains continuous data
%   cfg.channel                 = cell-array with channel labels, see FT_CHANNELSELECTION
%   cfg.comps                   = a vector with the components to plot (ex. 1:10) (optional)
%   cfg.zscale                  = [zmin zmax] or 'auto' (default = 'auto')
%   cfg.blocksize               = number (in seconds), only aplicable if data contains only 1 (long) trial
%   cfg.viewmode                = string, 'butterfly', 'vertical', 'component' (default = 'butterfly')
%   cfg.artfctdef.xxx.artifact  = Nx2 matrix with artifact segments see FT_ARTIFACT_xxx functions
%   cfg.selectfeature           = string, name of feature to be selected/added (default = 'visual')
%   cfg.selectmode              = string, what to do with a selection, can be 'mark', or 'eval' (default = 'mark')
%                                 'mark': artfctdef field is updated, 'eval': the function defined in
%                                 cfg.selfun is evaluated f.i. browse_movieplotER calls movieplotER which makes
%                                 a movie of the selected data
%   cfg.colorgroups             = 'sequential' 'labelcharx' (x = xth character in label), 'chantype' or
%                                  vector with length(data/hdr.label) defining groups (default = 'sequential')
%   cfg.channelcolormap         = COLORMAP (default = customized lines map with 15 colors)
%   cfg.selfun                  = string, name of function which is evaluated if selectmode is set to 'eval'.
%                                  The selected data and the selcfg are passed on to this function.
%   cfg.selcfg                  = configuration options for selfun
%   cfg.eegscale                = number, scaling to apply to the EEG channels prior to display
%   cfg.eogscale                = number, scaling to apply to the EOG channels prior to display
%   cfg.ecgscale                = number, scaling to apply to the ECG channels prior to display
%   cfg.emgscale                = number, scaling to apply to the EMG channels prior to display
%   cfg.megscale                = number, scaling to apply to the MEG channels prior to display
%
% The "artifact" field in the output cfg is a Nx2 matrix comparable to the
% "trl" matrix of FT_DEFINETRIAL. The first column of which specifying the
% beginsamples of an artifact period, the second column contains the
% endsamples of the artifactperiods.
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file.
% These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% NOTE for debugging: in case the databrowser crashes, use delete(gcf) to kill the figure.
%
% See also FT_PREPROCESSING, FT_REJECTARTIFACT, FT_ARTIFACT_EOG, FT_ARTIFACT_MUSCLE,
% FT_ARTIFACT_JUMP, FT_ARTIFACT_MANUAL, FT_ARTIFACT_THRESHOLD, FT_ARTIFACT_CLIP, FT_ARTIFACT_ECG

% Copyright (C) 2009, Robert Oostenveld, Ingrid Niewenhuis
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

% record start time and total processing time
ftFuncTimer = tic();
ftFuncClock = clock();
ftFuncMem   = memtic();

% set defaults for optional cfg.input and or cfg.outputfile
if ~isfield(cfg, 'inputfile'),       cfg.inputfile = [];               end
if ~isfield(cfg, 'outputfile'),      cfg.outputfile = [];              end

% set the defaults
if ~isfield(cfg, 'channel'),         cfg.channel = 'all';             end
if ~isfield(cfg, 'zscale'),          cfg.zscale = 'auto';             end
if ~isfield(cfg, 'artfctdef'),       cfg.artfctdef = struct;          end
if ~isfield(cfg, 'selectfeature'),   cfg.selectfeature = 'visual';    end % string or cell-array
if ~isfield(cfg, 'selectmode'),      cfg.selectmode = 'mark';         end
if ~isfield(cfg, 'viewmode'),        cfg.viewmode = 'butterfly';      end % butterfly, vertical, component, settings
if ~isfield(cfg, 'blocksize'),       cfg.blocksize = 1;               end % only for segmenting continuous data, i.e. one long trial
if ~isfield(cfg, 'preproc'),         cfg.preproc = [];                end % see preproc for options
if ~isfield(cfg, 'event'),           cfg.event = [];                  end
if ~isfield(cfg, 'selfun'),          cfg.selfun = 'browse_multiplotER';   end
if ~isfield(cfg, 'selcfg'),          cfg.selcfg = [];                     end
if ~isfield(cfg, 'colorgroups'),     cfg.colorgroups = 'sequential';      end
if ~isfield(cfg, 'channelcolormap'), cfg.channelcolormap = [0.75 0 0;0 0 1;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];   end
if ~isfield(cfg, 'eegscale'),        cfg.eegscale = [];                   end
if ~isfield(cfg, 'eogscale'),        cfg.eogscale = [];                   end
if ~isfield(cfg, 'ecgscale'),        cfg.ecgscale = [];                   end
if ~isfield(cfg, 'emgscale'),        cfg.emgscale = [];                   end
if ~isfield(cfg, 'megscale'),        cfg.megscale = [];                   end

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
  data = ft_checkdata(data, 'datatype', {'raw', 'comp'}, 'feedback', 'yes', 'hassampleinfo', 'yes');
  % fetch the header from memory
  hdr = ft_fetch_header(data);
  if ~isfield(cfg, 'continuous') && length(data.trial) == 1
    cfg.continuous = 'yes';
  elseif ~isfield(cfg, 'continuous') && length(data.trial) > 1
    cfg.continuous = 'no';
  end
else
  % check if the input cfg is valid for this function
  cfg = ft_checkconfig(cfg, 'dataset2files', {'yes'});
  cfg = ft_checkconfig(cfg, 'required', {'headerfile', 'datafile'});
  cfg = ft_checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
  cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});
  % read the header from file
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
  if ~isfield(cfg, 'continuous')
    if hdr.nTrials==1
      cfg.continuous = 'yes';
    else
      cfg.continuous = 'no';
    end
  end
end

if hasdata && isfield(data, 'topo') && strcmp(cfg.viewmode, 'component')
  if ~isfield(cfg, 'comp')
    cfg.comp = 1:10; % to avoid plotting 274 components topographically
  end
  cfg.channel = data.label(cfg.comp);
end

if ischar(cfg.selectfeature)
  % ensure that it is a cell array
  cfg.selectfeature = {cfg.selectfeature};
end

% get some initial parameters from the data
if hasdata
  
  % check whether data has been resampled
  if isfield(data, 'cfg') && isempty(ft_findcfg(data.cfg, 'origfs'))
    resampled = false;
  else
    resampled = true;
  end
  
  % fetch the events
  event = ft_fetch_event(data);
  
  cfg.channel = ft_channelselection(cfg.channel, data.label);
  chansel = match_str(data.label, cfg.channel);
  fsample = 1/mean(diff(data.time{1}));
  Nchans  = length(chansel);
  
  % this is how the input data is segmented
  trlorg = zeros(numel(data.trial), 3);
  trlorg(:,[1 2]) = data.sampleinfo;

  % recreate offset vector (databrowser depends on this for visualisation)
  for ntrl = 1:numel(data.trial)
    trlorg(ntrl,3) = time2offset(data.time{ntrl}, data.fsample);
  end
  Ntrials = size(trlorg, 1);
  
  if strcmp(cfg.viewmode, 'component')
    if ~isfield(cfg, 'layout')
      error('You need to specify a layout-file when browsing through components');
    end
    % read or create the layout that will be used for the topoplots
    cfg.layout = ft_prepare_layout(cfg, data);
    
    if ~isfield(cfg, 'comp')
      cfg.comp = 1:10; % to avoid plotting 274 components topographically
    end
    
    cfg.channel = data.label(cfg.comp);
  end
  
else
  % data has not been resampled
  resampled = false;
  
  % read the events
  if isempty(cfg.event)
    event = ft_read_event(cfg.dataset);
  else
    event = cfg.event;
  end
  
  cfg.channel = ft_channelselection(cfg.channel, hdr.label);
  chansel = match_str(hdr.label, cfg.channel);
  fsample = hdr.Fs;
  Nchans  = length(chansel);
  Ntrials = hdr.nTrials;
  
  % construct trl-matrix for data from file on disk
  trlorg = zeros(Ntrials,3);
  for k = 1:Ntrials
    trlorg(k,[1 2]) = [1 hdr.nSamples] + [hdr.nSamples hdr.nSamples] .* (k-1);
  end
  
end % if hasdata

if Nchans == 0
  error('no channels to display');
end

if Ntrials == 0
  error('no trials to display');
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
  chan_colors = [R(cfg.colorgroups(:)) G(cfg.colorgroups(:)) B(cfg.colorgroups(:))];
elseif strcmp(cfg.colorgroups, 'chantype')
  type = ft_chantype(labels_all);
  [tmp1 tmp2 cfg.colorgroups] = unique(type);
  fprintf('%3d colorgroups were identified\n',length(tmp1))
  R = cfg.channelcolormap(:,1);
  G = cfg.channelcolormap(:,2);
  B = cfg.channelcolormap(:,3);
  chan_colors = [R(cfg.colorgroups(:)) G(cfg.colorgroups(:)) B(cfg.colorgroups(:))];
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
  chan_colors = [R(cfg.colorgroups(:)) G(cfg.colorgroups(:)) B(cfg.colorgroups(:))];
elseif strcmp(cfg.colorgroups, 'sequential')
  % no grouping
  chan_colors = lines(length(labels_all));
else
  error('do not understand cfg.colorgroups')
end

% collect the artifacts that have been detected from cfg.artfctdef.xxx.artifact
artlabel = fieldnames(cfg.artfctdef);
sel      = zeros(size(artlabel));
artifact = cell(size(artlabel));
for i=1:length(artlabel)
  sel(i) = issubfield(cfg.artfctdef, [artlabel{i} '.artifact']);
  if sel(i)
    artifact{i} = getsubfield(cfg.artfctdef, [artlabel{i} '.artifact']);
    num         = size(artifact{i}, 1);
    if isempty(num)
      num = 0;
    end
    fprintf('detected %3d %s artifacts\n', num, artlabel{i});
  end
end
artifact = artifact(sel==1);
artlabel = artlabel(sel==1);

for i=1:length(cfg.selectfeature)
  if ~any(strcmp(cfg.selectfeature{i}, artlabel))
    artifact = {[],                   artifact{:}};
    artlabel = {cfg.selectfeature{i}, artlabel{:}};
  end
end

if length(artlabel) > 9
  error('only upto 9 artifacts groups supported')
end

% make artdata representing all artifacts in a "raw data" format
datendsample = max(trlorg(:,2));
artdat = convert_event(artifact, 'boolvec', 'endsample', datendsample);

artdata = [];
artdata.trial{1}       = artdat; % every artifact is a "channel"
artdata.time{1}        = offset2time(0, fsample, datendsample);
artdata.label          = artlabel;
artdata.fsample        = fsample;
artdata.cfg.trl        = [1 datendsample 0];

if ischar(cfg.zscale) && strcmp(cfg.zscale, 'auto')
  if nargin>1
    dat = data.trial{1}(chansel,:);
    time = data.time{1};
  else
    % one second of data is read from file to determine the vertical scaling
    dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', 1, 'endsample', round(hdr.Fs), 'chanindx', chansel, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat, 'headerformat', cfg.headerformat);
    time = (1:hdr.nSamples) / fsample;
  end
  
  minval = min(dat(:));
  maxval = max(dat(:));
  mintime = min(time(:));
  maxtime = max(time(:));
  cfg.zscale = max(abs(minval), abs(maxval));
  cfg.yscale = max(abs(mintime), abs(maxtime));
end

h = figure;
set(h, 'KeyPressFcn',           @keyboard_cb);
set(h, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonDownFcn'});
set(h, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonUpFcn'});
set(h, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonMotionFcn'});

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
  opt.orgdata   = [];      % this means that it will look in opt.cfg.dataset
end
opt.artdata  = artdata;
opt.cfg      = cfg;        % the configuration of this function, not of the preprocessing
opt.hdr      = hdr;
opt.event    = event;
opt.trlop    = 1;          % active trial being displayed
opt.ftsel    = find(strcmp(artlabel,cfg.selectfeature)); % current artifact/feature being selected
opt.trlorg   = trlorg;
opt.fsample  = fsample;
opt.artcol   = [0.9686 0.7608 0.7686; 0.7529 0.7098 0.9647; 0.7373 0.9725 0.6824;0.8118 0.8118 0.8118; 0.9725 0.6745 0.4784; 0.9765 0.9176 0.5686; 0.6863 1 1; 1 0.6863 1; 0 1 0.6000];
opt.chan_colors = chan_colors;
opt.cleanup  = false;      % this is needed for a corrent handling if the figure is closed (either in the corner or by "q")
opt.compindx = [];         % index of components to be drawn (if viewmode = "component")
opt.resampled = resampled;

if strcmp(cfg.continuous, 'yes')
  opt.trialname = 'segment';
else
  opt.trialname = 'trial';
end
guidata(h, opt);

% make the user interface elements for the data view
uicontrol('tag', 'group1', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', opt.trialname, 'userdata', 't')
uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'leftarrow')
uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'rightarrow')

uicontrol('tag', 'group1', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'channel','userdata', 'c')
uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'uparrow')
uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'downarrow')

uicontrol('tag', 'group1a', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'horizontal', 'userdata', 'h')
uicontrol('tag', 'group2a', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'shift+leftarrow')
uicontrol('tag', 'group2a', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'shift+rightarrow')
if strcmp(cfg.continuous, 'no')
  ft_uilayout(h, 'tag', 'group1a', 'visible', 'on', 'retag', 'group1');
  ft_uilayout(h, 'tag', 'group2a', 'visible', 'on', 'retag', 'group2');
else
  ft_uilayout(h, 'tag', 'group1a', 'visible', 'on', 'retag', 'group1');
  ft_uilayout(h, 'tag', 'group2a', 'visible', 'on', 'retag', 'group2');
end

uicontrol('tag', 'group1', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'vertical', 'userdata', 'v')
uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'shift+downarrow')
uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'shift+uparrow')

% legend artifacts/features
for iArt = 1:length(artlabel)
  uicontrol('tag', 'group3', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', artlabel{iArt}, 'userdata', num2str(iArt), 'position', [0.91, 0.9 - ((iArt-1)*0.09), 0.08, 0.04], 'backgroundcolor', opt.artcol(iArt,:))
  uicontrol('tag', 'group3', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', ['shift+' num2str(iArt)], 'position', [0.91, 0.85 - ((iArt-1)*0.09), 0.03, 0.04], 'backgroundcolor', opt.artcol(iArt,:))
  uicontrol('tag', 'group3', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', ['control+' num2str(iArt)], 'position', [0.96, 0.85 - ((iArt-1)*0.09), 0.03, 0.04], 'backgroundcolor', opt.artcol(iArt,:))
end

if strcmp(cfg.viewmode, 'butterfly')
  % button to find label of nearest channel to datapoint
  uicontrol('tag', 'group3', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'identify', 'userdata', 'i', 'position', [0.91, 0.1, 0.08, 0.05], 'backgroundcolor', [1 1 1])
end

ft_uilayout(h, 'tag', 'group1', 'width', 0.10, 'height', 0.05);
ft_uilayout(h, 'tag', 'group2', 'width', 0.05, 'height', 0.05);

ft_uilayout(h, 'tag', 'group1', 'style', 'pushbutton', 'callback', @keyboard_cb);
ft_uilayout(h, 'tag', 'group2', 'style', 'pushbutton', 'callback', @keyboard_cb);
ft_uilayout(h, 'tag', 'group3', 'style', 'pushbutton', 'callback', @keyboard_cb);

ft_uilayout(h, 'tag', 'group1', 'retag', 'viewui');
ft_uilayout(h, 'tag', 'group2', 'retag', 'viewui');
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
    opt = guidata(h);
    if opt.cleanup
      delete(h);
    end
  end
  
  % add the updated artifact definitions to the output cfg
  for i=1:length(opt.artdata.label)
    cfg.artfctdef.(opt.artdata.label{i}).artifact = convert_event(opt.artdata.trial{1}(i,:), 'artifact');
  end
end % if nargout

% add version information to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id = '$Id$';

% add information about the Matlab version used to the configuration
cfg.callinfo.matlab = version();
  
% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.procmem  = memtoc(ftFuncMem);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername();
fprintf('the call to "%s" took %d seconds and an estimated %d MB\n', mfilename, round(cfg.callinfo.proctime), round(cfg.callinfo.procmem/(1024*1024)));

% remember the configuration details of the input data
if hasdata && isfield(data, 'cfg')
  cfg.previous = data.cfg;
end

% remember the exact configuration details in the output
dataout.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', dataout); % use the variable name "data" in the output file
end

end % main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanup_cb(h, eventdata)
opt = guidata(h);
opt.cleanup = true;
guidata(h, opt);
uiresume
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function definetrial_cb(h, eventdata)
opt = guidata(h);
if strcmp(opt.cfg.continuous, 'no')
  % keep the original trial definition for visualisation
  opt.trlvis = opt.trlorg;
else
  % construct a trial definition for visualisation
  if isfield(opt, 'trlvis')
    thistrlbeg = opt.trlvis(opt.trlop,1);
    thistrlend = opt.trlvis(opt.trlop,2);
    % remember a representative sample of the current trial
    % thissample = round((thistrlbeg+thistrlend)/2);
    thissample = thistrlbeg;
  end
  % look at opt.cfg.blocksize and make opt.trl accordingly
  % if original data contains more than one trial, it will fail in ft_fetch_data
  datbegsample = min(opt.trlorg(:,1));
  datendsample = max(opt.trlorg(:,2));
  smppertrl  = round(opt.fsample * opt.cfg.blocksize);
  begsamples = datbegsample:smppertrl:datendsample;
  endsamples = datbegsample+smppertrl-1:smppertrl:datendsample;
  if numel(endsamples)<numel(begsamples)
    endsamples(end+1) = datendsample;
  end
  trlvis = [];
  trlvis(:,1) = begsamples';
  trlvis(:,2) = endsamples';
  if size(opt.trlorg,1) > 1 || isempty(opt.orgdata)
    % offset is now (re)defined that 1st sample is time 0
    trlvis(:,3) = begsamples-1;
  else
    % offset according to original time axis
    trlvis(:,3) = opt.trlorg(3) + begsamples - opt.trlorg(1);
  end
  
  if isfield(opt, 'trlvis')
    % update the current trial counter and try to keep the current sample the same
    % opt.trlop   = nearest(round((begsamples+endsamples)/2), thissample);
    opt.trlop   = nearest(begsamples, thissample);
  end
  opt.trlvis  = trlvis;
end % if continuous
guidata(h, opt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function help_cb(h, eventdata)
fprintf('------------------------------------------------------------------------------------\n')
fprintf('You can use the following buttons in the data viewer\n')
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
fprintf('q                  : quit\n');
fprintf('------------------------------------------------------------------------------------\n')
fprintf('\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_range_cb(range, h) %range 1X4 in sec relative to current trial
opt = guidata(h);

switch opt.cfg.viewmode
  case 'butterfly'
    begsample = opt.trlvis(opt.trlop,1);
    endsample = opt.trlvis(opt.trlop,2);
    offset    = opt.trlvis(opt.trlop,3);
    % determine the selection
    if strcmp(opt.trialname, 'trial')
      % this is appropriate when the offset is defined according to a
      % different trigger in each trial, which is usually the case in trial data
      begsel = round(range(1)*opt.fsample+begsample-offset-1);
      endsel = round(range(2)*opt.fsample+begsample-offset);
    elseif strcmp(opt.trialname, 'segment')
      % this is appropriate when the offset is defined according to a
      % one trigger, which is always the case in segment data [I think ingnie]
      begsel = round(range(1)*opt.fsample+1);
      endsel = round(range(2)*opt.fsample+1);
    end
    % the selection should always be confined to the current trial
    begsel = max(begsample, begsel);
    endsel = min(endsample, endsel);
    
  case {'vertical', 'component'}
    % the range should be in the displayed box
    range(1) = max(opt.hpos(1), range(1));
    range(2) = max(opt.hpos(1), range(2));
    range(1) = min(opt.hpos(2), range(1));
    range(2) = min(opt.hpos(2), range(2));
    range = (range - opt.hpos(1)) / (opt.hpos(2) - opt.hpos(1)); % left side of the box becomes 0, right side becomes 1
    range = range * (opt.hlim(2) - opt.hlim(1)) + opt.hlim(1);   % 0 becomes hlim(1), 1 becomes hlim(2)
    
    begsample = opt.trlvis(opt.trlop,1);
    endsample = opt.trlvis(opt.trlop,2);
    offset    = opt.trlvis(opt.trlop,3);
    % determine the selection
    if strcmp(opt.trialname, 'trial')
      % this is appropriate when the offset is defined according to a
      % different trigger in each trial, which is usually the case in trial data
      begsel = round(range(1)*opt.fsample+begsample-offset-1);
      endsel = round(range(2)*opt.fsample+begsample-offset);
    elseif strcmp(opt.trialname, 'segment')
      % this is appropriate when the offset is defined according to
      % one trigger, which is always the case in segment data [I think ingnie]
      begsel = round(range(1)*opt.fsample+1);
      endsel = round(range(2)*opt.fsample+1);
    end
    % the selection should always be confined to the current trial
    begsel = max(begsample, begsel);
    endsel = min(endsample, endsel);
    
  otherwise
    error('unknown opt.cfg.viewmode "%s"', opt.cfg.viewmode);
end % switch

if strcmp(opt.cfg.selectmode, 'disp')
  % FIXME this is only for debugging
  disp([begsel endsel])
  
elseif strcmp(opt.cfg.selectmode, 'mark')
  % mark or unmark artifacts
  artval = opt.artdata.trial{1}(opt.ftsel, begsel:endsel);
  artval = any(artval,1);
  if any(artval)
    fprintf('there is overlap with the active artifact (%s), disable this artifact\n',opt.artdata.label{opt.ftsel});
    opt.artdata.trial{1}(opt.ftsel, begsel:endsel) = 0;
  else
    fprintf('there is no overlap with the active artifact (%s), mark this as a new artifact\n',opt.artdata.label{opt.ftsel});
    opt.artdata.trial{1}(opt.ftsel, begsel:endsel) = 1;
  end
  
elseif strcmp(opt.cfg.selectmode, 'eval')
  % cut out the requested data segment
  seldata.label    = opt.curdat.label;
  seldata.time{1}  = offset2time(offset+begsel-begsample, opt.fsample, endsel-begsel+1);
  seldata.trial{1} = ft_fetch_data(opt.curdat, 'begsample', begsel, 'endsample', endsel);
  seldata.fsample  = opt.fsample;
  seldata.cfg.trl  = [begsel endsel offset];
  
  feval(opt.cfg.selfun, opt.cfg.selcfg, seldata);
  
else
  error('unknown value for cfg.selectmode');
end

guidata(h, opt);
redraw_cb(h);
end % function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keyboard_cb(h, eventdata)
opt = guidata(h);

if isempty(eventdata)
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
  h = getparent(h);
else
  % determine the key that was pressed on the keyboard
  key = parseKeyboardEvent(eventdata);
end

switch key
  case {'1' '2' '3' '4' '5' '6' '7' '8' '9'}
    % switch to another artifact type
    opt.ftsel = str2double(key);
    numart = size(opt.artdata.trial{1}, 1);
    if opt.ftsel > numart
      fprintf('data has no artifact type %i \n', opt.ftsel)
    else
      guidata(h, opt);
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
        guidata(h, opt);
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
        guidata(h, opt);
        redraw_cb(h, eventdata);
      end
    end
  case 'leftarrow'
    opt.trlop = max(opt.trlop - 1, 1); % should not be smaller than 1
    guidata(h, opt);
    redraw_cb(h, eventdata);
  case 'rightarrow'
    opt.trlop = min(opt.trlop + 1, size(opt.trlvis,1)); % should not be larger than the number of trials
    guidata(h, opt);
    redraw_cb(h, eventdata);
  case 'uparrow'
    chansel = match_str(opt.hdr.label, opt.cfg.channel);
    minchan = min(chansel);
    numchan = length(chansel);
    chansel = minchan - numchan : minchan - 1;
    if min(chansel)<1
      chansel = chansel - min(chansel) + 1;
    end
    % convert numeric array into cell-array with channel labels
    opt.cfg.channel = opt.hdr.label(chansel);
    disp(opt.cfg.channel);
    guidata(h, opt);
    redraw_cb(h, eventdata);
  case 'downarrow'
    chansel = match_str(opt.hdr.label, opt.cfg.channel);
    maxchan = max(chansel);
    numchan = length(chansel);
    chansel = maxchan + 1 : maxchan + numchan;
    if max(chansel)>length(opt.hdr.label)
      chansel = chansel - (max(chansel) - length(opt.hdr.label));
    end
    % convert numeric array into cell-array with channel labels
    opt.cfg.channel = opt.hdr.label(chansel);
    disp(opt.cfg.channel);
    guidata(h, opt);
    redraw_cb(h, eventdata);
  case 'shift+leftarrow'
    opt.cfg.blocksize = opt.cfg.blocksize*sqrt(2);
    disp(opt.cfg.blocksize);
    guidata(h, opt);
    definetrial_cb(h, eventdata);
    redraw_cb(h, eventdata);
  case 'shift+rightarrow'
    opt.cfg.blocksize = opt.cfg.blocksize/sqrt(2);
    disp(opt.cfg.blocksize);
    guidata(h, opt);
    definetrial_cb(h, eventdata);
    redraw_cb(h, eventdata);
  case 'shift+uparrow'
    opt.cfg.zscale = opt.cfg.zscale/sqrt(2);
    guidata(h, opt);
    redraw_cb(h, eventdata);
  case 'shift+downarrow'
    opt.cfg.zscale = opt.cfg.zscale*sqrt(2);
    guidata(h, opt);
    redraw_cb(h, eventdata);
  case 'q'
    guidata(h, opt);
    cleanup_cb(h);
  case 't'
    % select the trial to display
    response = inputdlg(sprintf('%s to display', opt.trialname), 'specify', 1, {num2str(opt.trlop)});
    if ~isempty(response)
      opt.trlop = str2double(response);
      opt.trlop = min(opt.trlop, size(opt.trlvis,1)); % should not be larger than the number of trials
      opt.trlop = max(opt.trlop, 1); % should not be smaller than 1
    end
    guidata(h, opt);
    redraw_cb(h, eventdata);
  case 'h'
    % select the horizontal scaling
    response = inputdlg('horizontal scale', 'specify', 1, {num2str(opt.cfg.blocksize)});
    if ~isempty(response)
      opt.cfg.blocksize = str2double(response);
    end
    guidata(h, opt);
    definetrial_cb(h, eventdata);
    redraw_cb(h, eventdata);
  case 'v'
    % select the vertical scaling
    response = inputdlg('vertical scale, number or ''maxabs'')', 'specify', 1, {num2str(opt.cfg.zscale)});
    if ~isempty(response)
      if isnan(str2double(response)) && strcmp(response, 'maxabs')
        minval = min(opt.curdat.trial{1}(:));
        maxval = max(opt.curdat.trial{1}(:));
        opt.cfg.zscale = max(abs([minval maxval]));
      else
        opt.cfg.zscale = str2double(response);
      end
    end
    guidata(h, opt);
    redraw_cb(h, eventdata);
  case 'c'
    % select channels
    select = match_str(opt.hdr.label, opt.cfg.channel);
    select = select_channel_list(opt.hdr.label, select);
    opt.cfg.channel = opt.hdr.label(select);
    guidata(h, opt);
    redraw_cb(h, eventdata);
  case 'i'
    if strcmp(opt.cfg.viewmode, 'butterfly')
      % click in data and get name of nearest channel
      fprintf('click in the figure to identify the name of the closest channel\n');
      val = ginput(1);
      channame = val2nearestchan(opt.curdat,val);
      channb = match_str(opt.curdat.label,channame);
      fprintf('channel name: %s\n',channame);
      redraw_cb(h, eventdata);
      hold on
      xtext = opt.cfg.zscale - 0.1*opt.cfg.zscale;
      ft_plot_text(val(1), xtext, channame, 'FontSize', 16);
      plot(opt.curdat.time{1}, opt.curdat.trial{1}(channb,:),'k','LineWidth',2)
    end
  case 'control+control'
    % do nothing
  case 'shift+shift'
    % do nothing
  case 'alt+alt'
    % do nothing
  otherwise
    guidata(h, opt);
    help_cb(h);
end
uiresume(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toggle_viewmode_cb(h, eventdata, varargin)
opt = guidata(getparent(h));
if ~isempty(varargin) && ischar(varargin{1})
  opt.cfg.viewmode = varargin{1};
end
guidata(getparent(h), opt);
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
opt = guidata(h);
figure(h); % ensure that the calling figure is in the front

fprintf('redrawing with viewmode %s\n', opt.cfg.viewmode);

begsample = opt.trlvis(opt.trlop, 1);
endsample = opt.trlvis(opt.trlop, 2);
offset    = opt.trlvis(opt.trlop, 3);
chanindx  = match_str(opt.hdr.label, opt.cfg.channel);

if ~isempty(opt.event)
  % select only the events in the current time window
  event     = opt.event;
  evtsample = [event(:).sample];
  event     = event(evtsample>=begsample & evtsample<=endsample);
else
  event = [];
end

if isempty(opt.orgdata)
  fprintf('reading data... ');
  dat = ft_read_data(opt.cfg.datafile, 'header', opt.hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', strcmp(opt.cfg.continuous, 'no'), 'dataformat', opt.cfg.dataformat, 'headerformat', opt.cfg.headerformat);
else
  fprintf('fetching data... ');
  dat = ft_fetch_data(opt.orgdata, 'header', opt.hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx);
end
fprintf('done\n');

fprintf('fetching artifacts... ');
art = ft_fetch_data(opt.artdata, 'begsample', begsample, 'endsample', endsample);
fprintf('done\n');

% apply preprocessing and determine the time axis
fprintf('preprocessing data... ');
[dat, lab, tim] = preproc(dat, opt.hdr.label(chanindx), offset2time(offset, opt.fsample, size(dat,2)), opt.cfg.preproc);
fprintf('done\n');

opt.curdat.label    = lab;
opt.curdat.time{1}  = tim;
opt.curdat.trial{1} = dat;
opt.curdat.fsample  = opt.fsample;
opt.curdat.cfg.trl  = [begsample endsample offset];

% apply scaling to selected channels 
% using wildcard to support subselection of channels
if ~isempty(opt.cfg.megscale)
  chansel = match_str(lab, ft_channelselection('MEG*', lab));
  dat(chansel,:) = dat(chansel,:) .* opt.cfg.megscale;
end
if ~isempty(opt.cfg.eegscale)
  chansel = match_str(lab, ft_channelselection('EEG*', lab));
  dat(chansel,:) = dat(chansel,:) .* opt.cfg.eegscale;
end
if ~isempty(opt.cfg.eogscale)
  chansel = match_str(lab, ft_channelselection('EOG*', lab));
  dat(chansel,:) = dat(chansel,:) .* opt.cfg.eogscale;
end
if ~isempty(opt.cfg.ecgscale)
  chansel = match_str(lab, ft_channelselection('ECG*', lab));
  dat(chansel,:) = dat(chansel,:) .* opt.cfg.ecgscale;
end
if ~isempty(opt.cfg.emgscale)
  chansel = match_str(lab, ft_channelselection('EMG*', lab));
  dat(chansel,:) = dat(chansel,:) .* opt.cfg.emgscale;
end

fprintf('plotting data... ');
switch opt.cfg.viewmode
  case 'butterfly'
    cla;       % clear the content in the current axis
    % to assure current feature is plotted on top
    ordervec = 1:length(opt.artdata.label);
    ordervec(opt.ftsel) = [];
    ordervec(end+1) = opt.ftsel;
    
    for i = ordervec
      tmp = diff([0 art(i,:) 0]);
      artbeg = find(tmp==+1);
      artend = find(tmp==-1) - 1;
      for j=1:numel(artbeg)
        ft_plot_box([tim(artbeg(j)) tim(artend(j)) -opt.cfg.zscale opt.cfg.zscale], 'facecolor', opt.artcol(i,:), 'edgecolor', 'none');
      end
    end
    
    h_event = zeros(0, length(event));
    h_event_txt = zeros(0, length(event));
    if ~opt.resampled
      try
        % plot a line with text for each event
        for k=1:length(event)
          try
            eventstr = sprintf('%s=%s', event(k).type, num2str(event(k).value)); %value can be both number and string
          catch
            eventstr = 'unknown';
          end
          eventtim = (event(k).sample-begsample+offset)/opt.fsample;
          eventtim = (eventtim - opt.hlim(1)) / (opt.hlim(2) - opt.hlim(1));   % convert to value relative to box, i.e. from 0 to 1
          eventtim = eventtim * (opt.hpos(2) - opt.hpos(1)) + opt.hpos(1);     % convert from relative to actual value along the horizontal figure axis
          h_event(k) = ft_plot_line([eventtim eventtim], [0 1]);
          %       h_event(k) = ft_plot_line([eventtim eventtim], [-opt.cfg.zscale opt.cfg.zscale]);
          h_event_txt(k) = ft_plot_text(eventtim, ax(4)-0.01, eventstr);
        end
      end % try
    else
      if isfield(opt, 'orgdata') && isfield(opt.orgdata, 'sampleinfo') && isfield(opt.orgdata, 'offset')
        % find trials within this segment
        trlindx = find(((opt.orgdata.sampleinfo(:, 1)-opt.orgdata.offset) >= begsample & (opt.orgdata.sampleinfo(:, 1)-opt.orgdata.offset) <= endsample)==1);
        for t = 1:numel(trlindx)
          eventtim = (opt.orgdata.sampleinfo(trlindx(t), 1)-opt.orgdata.offset(trlindx(t))-begsample+offset)/opt.fsample;
          eventtim = (eventtim - opt.hlim(1)) / (opt.hlim(2) - opt.hlim(1));   % convert to value relative to box, i.e. from 0 to 1
          eventtim = eventtim * (opt.hpos(2) - opt.hpos(1)) + opt.hpos(1);     % convert from relative to actual value along the horizontal figure axis
          h_event(end+1) = ft_plot_line([eventtim eventtim], [-opt.cfg.zscale 1]);
          h_event_txt(end+1) = ft_plot_text(eventtim, ax(4)+.01, 'stim');
        end
      end
    end
        
    set(h_event, 'tag', 'events');
    set(h_event_txt, 'tag', 'events');
    set(gca,'ColorOrder',opt.chan_colors(chanindx,:)) % plot vector does not clear axis, therefore this is possible
    
    % plot the data on top of the box
    h_act = ft_plot_vector(tim, dat);
    set(h_act, 'tag', 'activations');
    
    ax(1) = tim(1);
    ax(2) = tim(end);
    ax(3) = -opt.cfg.zscale;
    ax(4) =  opt.cfg.zscale;
    axis(ax);
    title(sprintf('%s %d, time from %g to %g s', opt.trialname, opt.trlop, tim(1), tim(end)));
    xlabel('time');
    
    % set tags
    
    
  case 'vertical'
    cla;       % clear the content in the current axis
    tmpcfg = [];
    tmpcfg.layout = 'vertical';
    tmpcfg.channel = opt.cfg.channel;
    tmpcfg.skipcomnt = 'yes';
    tmpcfg.skipscale = 'yes';
    laytime = ft_prepare_layout(tmpcfg, opt.orgdata);
    
    hlim = [tim(1) tim(end)];
    vlim = [-opt.cfg.zscale +opt.cfg.zscale];
    
    % determine the position of each of the labels
    labelx = laytime.pos(:,1) - laytime.width/2 - 0.01;
    labely = laytime.pos(:,2);
    
    ax(1) = min(laytime.pos(:,1) - laytime.width/2);
    ax(2) = max(laytime.pos(:,1) + laytime.width/2);
    ax(3) = min(laytime.pos(:,2) - laytime.height/2);
    ax(4) = max(laytime.pos(:,2) + laytime.height/2);
    axis(ax)
    
    % remember the scaling of the horizontal axis, this is needed for mouse input
    hpos(1) = laytime.pos(1,1) - laytime.width(1)/2; % the position of the left  side of the timecourse box
    hpos(2) = laytime.pos(1,1) + laytime.width(1)/2; % the position of the right side of the timecourse box
    opt.hlim = hlim;
    opt.hpos = hpos;
    
    % to assure current feature is plotted on top
    ordervec = 1:length(opt.artdata.label);
    ordervec(opt.ftsel) = [];
    ordervec(end+1) = opt.ftsel;
    
    for j = ordervec
      tmp = diff([0 art(j,:) 0]);
      artbeg = find(tmp==+1);
      artend = find(tmp==-1) - 1;
      arttim = [tim(artbeg)' tim(artend)'];                            % convert the artifact sample number to time
      arttim = (arttim - opt.hlim(1)) / (opt.hlim(2) - opt.hlim(1));   % convert to value relative to box, i.e. from 0 to 1
      arttim = arttim * (opt.hpos(2) - opt.hpos(1)) + opt.hpos(1);     % convert from relative to actual value along the horizontal figure axis
      
      for k=1:numel(artbeg)
        ft_plot_box([arttim(k,1) arttim(k,2) ax(3) ax(4)], 'facecolor', opt.artcol(j,:), 'edgecolor', 'none');
      end
    end % for each of the artifact channels
    
    h_event = zeros(0, length(event));
    h_event_txt = zeros(0, length(event));
    if ~opt.resampled
      % plot a line with text for each event
      for k=1:length(event)
        try
          eventstr = sprintf('%s=%s', event(k).type, num2str(event(k).value)); %value can be both number and string
        catch
          eventstr = 'unknown';
        end
        eventtim = (event(k).sample-begsample+offset)/opt.fsample;
        eventtim = (eventtim - opt.hlim(1)) / (opt.hlim(2) - opt.hlim(1));   % convert to value relative to box, i.e. from 0 to 1
        eventtim = eventtim * (opt.hpos(2) - opt.hpos(1)) + opt.hpos(1);     % convert from relative to actual value along the horizontal figure axis
        h_event(k) = ft_plot_line([eventtim eventtim], [0 1]);
%       h_event(k) = ft_plot_line([eventtim eventtim], [-opt.cfg.zscale opt.cfg.zscale]);
        h_event_txt(k) = ft_plot_text(eventtim, ax(4)-0.01, eventstr);
      end
    else
      if isfield(opt, 'orgdata') && isfield(opt.orgdata, 'sampleinfo') && isfield(opt.orgdata, 'offset')
        % find trials within this segment
        trlindx = find(((opt.orgdata.sampleinfo(:, 1)-opt.orgdata.offset) >= begsample & (opt.orgdata.sampleinfo(:, 1)-opt.orgdata.offset) <= endsample)==1);
        for t = 1:numel(trlindx)
          eventtim = (opt.orgdata.sampleinfo(trlindx(t), 1)-opt.orgdata.offset(trlindx(t))-begsample+offset)/opt.fsample;
          eventtim = (eventtim - opt.hlim(1)) / (opt.hlim(2) - opt.hlim(1));   % convert to value relative to box, i.e. from 0 to 1
          eventtim = eventtim * (opt.hpos(2) - opt.hpos(1)) + opt.hpos(1);     % convert from relative to actual value along the horizontal figure axis
          h_event(end+1) = ft_plot_line([eventtim eventtim], [-opt.cfg.zscale 1]);
          h_event_txt(end+1) = ft_plot_text(eventtim, ax(4)+.01, 'stim');
        end
      end
    end
    % set tags
    set(h_event, 'tag', 'events');
    set(h_event_txt, 'tag', 'events');
    
    for i = 1:length(chanindx)
      datsel = i;
      laysel = match_str(laytime.label, opt.hdr.label(chanindx(i)));
      if ~isempty(datsel)
        h_text = ft_plot_text(labelx(laysel), labely(laysel), opt.hdr.label(chanindx(i)), 'HorizontalAlignment', 'right');
        h_act = ft_plot_vector(tim, dat(datsel, :), 'hpos', laytime.pos(laysel,1), 'vpos', laytime.pos(laysel,2), 'width', laytime.width(laysel), 'height', laytime.height(laysel), 'hlim', hlim, 'vlim', vlim, 'box', false, 'color', opt.chan_colors(chanindx(i),:));
      end
    end
    % set tags
    set(h_text, 'tag', 'activations');
    set(h_act, 'tag', 'activations');
    
    nticks = 11;
    set(gca, 'xTick', linspace(ax(1), ax(2), nticks))
    xTickLabel = cellstr(num2str( linspace(tim(1), tim(end), nticks)' , '%1.2f'))';
    set(gca, 'xTickLabel', xTickLabel)
    
    set(gca, 'yTick', [])
    if length(chanindx)<7
      % two ticks per channel
      set(gca, 'yTick', sort([laytime.pos(:,2)+(laytime.height(laysel)/2); laytime.pos(:,2)+(laytime.height(laysel)/4); laytime.pos(:,2)-(laytime.height(laysel)/4); laytime.pos(:,2)-(laytime.height(laysel)/2)]))
      yTickLabel = {num2str(-vlim(2)), num2str(-vlim(2)/2), num2str(vlim(2)/2), num2str(vlim(2))};
    elseif length(chanindx)> 6 && length(chanindx)< 20
      % one tick per channel
      set(gca, 'yTick', sort([laytime.pos(:,2)+(laytime.height(laysel)/4); laytime.pos(:,2)-(laytime.height(laysel)/4)]))
      yTickLabel = {num2str(-vlim(2)/2), num2str(vlim(2)/2)};
    else
      % no space for xticks
      yTickLabel = [];
    end
    tmp = yTickLabel;
    for chanloop = 2:length(chanindx)
      yTickLabel = [yTickLabel tmp];
    end
    set(gca, 'yTickLabel', yTickLabel)
    
    title(sprintf('%s %d, time from %g to %g s', opt.trialname, opt.trlop, tim(1), tim(end)));
    
  case 'component'
    % delete time courses
    delete(findobj(h,'tag', 'activations'));    
    delete(findobj(h,'tag', 'events'));
    delete(findobj(h,'tag', 'artifacts'));
    compindx = chanindx;
    clear chanindx
    
    tmpcfg = [];
    tmpcfg.layout = 'vertical';
    tmpcfg.channel = opt.cfg.channel;
    tmpcfg.skipcomnt = 'yes';
    tmpcfg.skipscale = 'yes';
    laytime = ft_prepare_layout(tmpcfg, opt.orgdata);
    
    ax(1) = min(laytime.pos(:,1) - laytime.width/2);
    ax(2) = max(laytime.pos(:,1) + laytime.width/2);
    ax(3) = min(laytime.pos(:,2) - laytime.height/2);
    ax(4) = max(laytime.pos(:,2) + laytime.height/2);
    
    tmpcfg = [];
    tmpcfg.layout = opt.cfg.layout;
    laychan = ft_prepare_layout(tmpcfg, opt.orgdata);
    
    
    
    % determine the position of each of the topographies
    laytopo.pos(:,1)  = laytime.pos(:,1) - laytime.width/2 - laytime.height*2;
    laytopo.pos(:,2)  = laytime.pos(:,2);
    laytopo.width     = laytime.height;
    laytopo.height    = laytime.height;
    laytopo.label     = laytime.label;
    
    % determine the position of each of the labels
    labelx = laytopo.pos(:,1) + laytopo.width;
    labely = laytopo.pos(:,2);
    
    hlim = [tim(1) tim(end)];
    vlim = [-opt.cfg.zscale +opt.cfg.zscale];
    hpos(1) = laytime.pos(1,1) - laytime.width(1)/2; % the position of the left  side of the timecourse box
    hpos(2) = laytime.pos(1,1) + laytime.width(1)/2; % the position of the right side of the timecourse box
    opt.hlim = hlim;
    opt.hpos = hpos;
    
        
    % check if topographies need to be redrawn
    redraw_topo = false;
    if ~isequal(opt.compindx, compindx)
      redraw_topo = true;
      cla; % clear axis
    end
    
    
    % to assure current feature is plotted on top
    ordervec = 1:length(opt.artdata.label);
    ordervec(opt.ftsel) = [];
    ordervec(end+1) = opt.ftsel;
    
    h_art = cell(1, ordervec);
    for j = ordervec
      tmp = diff([0 art(j,:) 0]);
      artbeg = find(tmp==+1);
      artend = find(tmp==-1) - 1;
      arttim = [tim(artbeg)' tim(artend)'];                            % convert the artifact sample number to time
      arttim = (arttim - opt.hlim(1)) / (opt.hlim(2) - opt.hlim(1));   % convert to value relative to box, i.e. from 0 to 1
      arttim = arttim * (opt.hpos(2) - opt.hpos(1)) + opt.hpos(1);     % convert from relative to actual value along the horizontal figure axis
      h_art{j} = zeros(1, length(artbeg));
      for k=1:numel(artbeg)
        h_art{j}(k) = ft_plot_box([arttim(k,1) arttim(k,2) ax(3) ax(4)], 'facecolor', opt.artcol(j,:), 'edgecolor', 'none');
      end
    end % for each of the artifact channels
    
    
    [sel1, sel2] = match_str(opt.orgdata.topolabel, laychan.label);
    chanx = laychan.pos(sel2,1);
    chany = laychan.pos(sel2,2);
    


    
    h_act = zeros(1, length(compindx));
    h_text = zeros(1, length(compindx));
    for i=1:length(compindx)
      datsel = i;
      laysel = match_str(laytime.label,opt.hdr.label(compindx(i)));
      if ~isempty(datsel)
        % plot the timecourse of this component
        h_act(i) = ft_plot_vector(tim, dat(datsel, :), 'hpos', laytime.pos(laysel,1), 'vpos', laytime.pos(laysel,2), 'width', laytime.width(laysel), 'height', laytime.height(laysel), 'hlim', hlim, 'vlim', vlim);
        
        if redraw_topo
          h_text(i) = ft_plot_text(labelx(laysel), labely(laysel), opt.hdr.label(compindx(i)));
          
          % plot the topography of this component
          chanz = opt.orgdata.topo(sel1,compindx(i));
          ft_plot_topo(chanx, chany, chanz./max(abs(chanz)), 'hpos', laytopo.pos(laysel,1), ...
            'vpos', laytopo.pos(laysel,2), 'mask', laychan.mask, ...
            'interplim', 'mask', 'outline', laychan.outline,  ...
            'width', laytopo.width(laysel), 'height', laytopo.height(laysel));
        end
        
        axis equal
        drawnow
      end
    end
    
    % draw events    
    h_event = zeros(0, length(event));
    h_event_txt = zeros(0, length(event));
    if ~opt.resampled
      % plot a line with text for each event
      for k=1:length(event)
        try
          eventstr = sprintf('%s=%s', event(k).type, num2str(event(k).value)); %value can be both number and string
        catch
          eventstr = 'unknown';
        end
        eventtim = (event(k).sample-begsample+offset)/opt.fsample;
        eventtim = (eventtim - opt.hlim(1)) / (opt.hlim(2) - opt.hlim(1));   % convert to value relative to box, i.e. from 0 to 1
        eventtim = eventtim * (opt.hpos(2) - opt.hpos(1)) + opt.hpos(1);     % convert from relative to actual value along the horizontal figure axis
        h_event(k) = ft_plot_line([eventtim eventtim], [0 1]);
%       h_event(k) = ft_plot_line([eventtim eventtim], [-opt.cfg.zscale opt.cfg.zscale]);
        h_event_txt(k) = ft_plot_text(eventtim, ax(4)-0.01, eventstr);
      end
    else
      if isfield(opt, 'orgdata') && isfield(opt.orgdata, 'sampleinfo') && isfield(opt.orgdata, 'offset')
        % find trials within this segment
        trlindx = find(((opt.orgdata.sampleinfo(:, 1)-opt.orgdata.offset) >= begsample & (opt.orgdata.sampleinfo(:, 1)-opt.orgdata.offset) <= endsample)==1);
        for t = 1:numel(trlindx)
          eventtim = (opt.orgdata.sampleinfo(trlindx(t), 1)-opt.orgdata.offset(trlindx(t))-begsample+offset)/opt.fsample;
          eventtim = (eventtim - opt.hlim(1)) / (opt.hlim(2) - opt.hlim(1));   % convert to value relative to box, i.e. from 0 to 1
          eventtim = eventtim * (opt.hpos(2) - opt.hpos(1)) + opt.hpos(1);     % convert from relative to actual value along the horizontal figure axis
          h_event(end+1) = ft_plot_line([eventtim eventtim], [-opt.cfg.zscale 1]);
          h_event_txt(end+1) = ft_plot_text(eventtim, ax(4)+.01, 'stim');
        end
      end
    end
    % set tags
    set(h_event, 'tag', 'events');
    set(h_event_txt, 'tag', 'events');
    
    set(h_act, 'tag', 'activations');
    for k = 1:numel(h_art)
        set(h_art{k}, 'tag', 'artifacts');
    end
    
    h_topo = findobj(h, 'type', 'surface');
    set(h_text, 'tag', 'comptopo')
    set(h_topo, 'tag', 'comptopo')
    
    opt.compindx = compindx;
    
    set(gca, 'xTick', [])
    set(gca, 'yTick', [])
    title(sprintf('%s %d, time from %g to %g s', opt.trialname, opt.trlop, tim(1), tim(end)));
    
    ax(1) = min(laytopo.pos(:,1) - laytopo.width/2);
    ax(2) = max(laytime.pos(:,1) + laytime.width/2);
    ax(3) = min(laytime.pos(:,2) - laytime.height/2);
    ax(4) = max(laytime.pos(:,2) + laytime.height/2);
    axis(ax)
    
    % remember the scaling of the horizontal axis, this is needed for mouse input
    hpos(1) = laytime.pos(1,1) - laytime.width(1)/2; % the position of the left  side of the timecourse box
    hpos(2) = laytime.pos(1,1) + laytime.width(1)/2; % the position of the right side of the timecourse box
    opt.hlim = hlim;
    opt.hpos = hpos;
    
  case 'settings'
    % FIXME implement further details
    
  otherwise
    error('unknown viewmode "%s"', opt.cfg.viewmode);
end % switch viewmode
fprintf('done\n');

guidata(h, opt);
end

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
  
end

